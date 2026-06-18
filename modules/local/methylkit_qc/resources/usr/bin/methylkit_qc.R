#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(methylKit)
  library(genomation)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})


option_list = list(
  make_option("--cov_dir", type="character", help="Directory with Bismark .cov files"),
  make_option("--metadata", type="character", help="Metadata CSV/TSV file"),
  make_option("--sample_suffix", type="character",
              default="_bismark.cov.gz",
              help="Suffix to strip from coverage filenames"),
  make_option("--group_column", type="character",
              help="Metadata column defining groups"),
  make_option("--group_case", type="character",
              help="Value defining treatment group (coded as 1)"),
  make_option("--assembly", type="character", default="mm10"),
  make_option("--cores", type="integer", default=4),
  make_option("--lo_count", type="integer", default=5),
  make_option("--lo_perc", type="double", default=NA),
  make_option("--hi_count", type="integer", default=NA),
  make_option("--hi_perc", type="double", default=NA),
  make_option("--destrand", type="logical", default=FALSE),
  make_option("--min_per_group", type = "character", default = NA)
)

opt = parse_args(OptionParser(option_list=option_list))

cov_dir <- opt$cov_dir        
metadata <- opt$metadata      
sample_suffix <- opt$sample_suffix 
group_column  <- opt$group_column   
group_case    <- opt$group_case    
assembly <- opt$assembly      
cores    <- opt$cores      

lo_count      <- opt$lo_count
lo_perc       <- if (is.na(opt$lo_perc)) NULL else opt$lo_perc
hi_count      <- if (is.na(opt$hi_count)) NULL else opt$hi_count
hi_perc       <- if (is.na(opt$hi_perc)) NULL else opt$hi_perc
destrand      <- opt$destrand
min_per_group <- opt$min_per_group

if (is.na(min_per_group) || min_per_group %in% c("NA", "NULL", "null", "")) {
  min_per_group <- NULL
} else {
  min_per_group <- as.integer(min_per_group)
  if (is.na(min_per_group)) {
    stop("--min_per_group must be an integer or NA")
  }
}

#input checks
stopifnot(
  dir.exists(cov_dir),
  file.exists(metadata),
  !is.null(group_column),
  !is.null(group_case)
)

#input

files_list = list.files(
  cov_dir,
  pattern = "\\.cov",
  full.names = TRUE
)

meta = read.delim(metadata, row.names = 1, check.names = FALSE, sep = ",")

sample_ids = sub(sample_suffix, "", basename(files_list))
files_list = files_list[sample_ids %in% rownames(meta)]
sample_ids = sample_ids[sample_ids %in% rownames(meta)]

if (length(files_list) == 0) {
  stop("No coverage files matched metadata row names after stripping sample_suffix")
}

meta = meta[sample_ids, , drop = FALSE]
files_list = as.list(files_list)

if (!group_column %in% colnames(meta)) {
  stop("group_column not found in metadata: ", group_column)
}

treat = ifelse(meta[[group_column]] == group_case, 1, 0)



# data normalization

methData = methRead(location=files_list,
                    sample.id=as.list(rownames(meta)), 
                    assembly=assembly, 
                    pipeline="bismarkCoverage", 
                    header=F, 
                    context="CpG", 
                    resolution="base", 
                    treatment=treat, 
                    mincov=5, 
                    dbtype="tabix", 
                    dbdir="methylDB")


methData.filt <- filterByCoverage(
  methData,
  lo.count = lo_count,
  lo.perc  = lo_perc,
  hi.count = hi_count,
  hi.perc  = hi_perc
)


methData.norm = normalizeCoverage(methData.filt, method = "median")

#methData.unite = unite(
#  methData.norm,
#  destrand = destrand,
#  mc.cores = cores,
#  min.per.group = min_per_group
#)

unite_args <- list(
  object = methData.norm,
  destrand = destrand,
  mc.cores = cores
)

if (!is.null(min_per_group)) {
  unite_args$min.per.group <- min_per_group
}

methData.unite <- do.call(unite, unite_args)

if (nrow(getData(methData.unite)) == 0) {
  stop("unite() returned zero CpG sites")
}

save(methData, methData.filt, methData.norm, methData.unite, meta, file = "meth_merged_data.rda")

# QC 

pdf("coverageStats_histogram.pdf")
for (i in seq_along(methData)) {
  getCoverageStats(methData[[i]], plot = TRUE)
}
dev.off()

pdf("methylationStats_histogram.pdf")
for (i in seq_along(methData)) {
  getMethylationStats(methData[[i]], plot = TRUE)
}
dev.off()

# PCA
plot_pca = function(df, pcx, pcy, color_var, var_exp) {
  ggplot(df, aes(.data[[pcx]], .data[[pcy]], color=.data[[color_var]])) +
    geom_point(size=3) +
    geom_text_repel(aes(label=sample), size=3) +
    theme_light() +
    labs(
      x = paste0(pcx, " (", var_exp[as.numeric(sub("PC","",pcx))], "%)"),
      y = paste0(pcy, " (", var_exp[as.numeric(sub("PC","",pcy))], "%)")
    )
}

my_prcomp = PCASamples(methData.unite, obj.return = TRUE)

var_exp = round(
  (my_prcomp$sdev^2) / sum(my_prcomp$sdev^2) * 100,
  1
)

var_df = data.frame(
  PC = factor(paste0("PC", seq_along(var_exp)),
              levels = paste0("PC", seq_along(var_exp))),
  variance = var_exp
)

scores = as.data.frame(my_prcomp$x)
#pca_df = cbind(scores[,1:4], meta) # TODO: why from 1 to 4?
pca_df = cbind(scores, meta)
pca_df$sample = rownames(pca_df)

pdf("PCA_covariates.pdf")

ggplot(var_df, aes(PC, variance)) +
  geom_col() +
  geom_text(aes(label=paste0(variance, "%")), vjust=-0.3) +
  theme_light()

for (col in colnames(meta)) {
  print(plot_pca(pca_df, "PC1", "PC2", col, var_exp))
}
#for (col in colnames(meta)) {
#  print(plot_pca(pca_df, "PC3", "PC4", col, var_exp))
#}

if (all(c("PC3", "PC4") %in% colnames(pca_df))) {
  for (col in colnames(meta)) {
    print(plot_pca(pca_df, "PC3", "PC4", col, var_exp))
  }
}

dev.off()

#DEV

#Plot common methylated sites
covered_sites <- sapply(methData.norm, function(x) nrow(getData(x)))
sample_names  <- sapply(methData.norm, function(x) x@sample.id)

cov_df <- data.frame(
  sample = sample_names,
  covered_sites = covered_sites,
  stringsAsFactors = FALSE
)

common_sites <- nrow(getData(methData.unite))

library(scales)

min_per_group_label <- if (is.null(min_per_group)) "All samples" else as.character(min_per_group)

p <- ggplot(cov_df, aes(x = sample, y = covered_sites)) +
  geom_col(fill = "steelblue") +
  geom_hline(
    yintercept = common_sites,
    color = "red",
    linewidth = 1
  ) +
  scale_y_continuous(labels = label_comma()) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = paste0("CpG site coverage per sample"),
    subtitle = paste0("Horizontal line shows common CpG sites used for PCA (", common_sites, ") - Replicate requirement : ", min_per_group_label,),
    x = "Sample",
    y = "Number of CpG sites"
  )

# Aggiungere altro plot oltre a quello che c'è già per mostrare quante posizioni hanno coverage>10
# in ciascuna delle posizioni dello unite, per ogni sample
# Così si capisce se nella definizione dello unite il contributo dei diversi samples è stato bilanciato o no

ggsave("covered_sites.pdf", p, width=8, height=5)
dev.off()

###  Clustering (Extra)
pdf("clustering_dendrogram.pdf")
par(mar=c(7,3,3,3))

mat <- methylKit::percMethylation(methData.unite)
corr <- cor(mat, use = "pairwise.complete.obs")

library(pheatmap)

pheatmap(
  corr,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10
)

clusterSamples(methData.unite, dist="correlation", method="ward", plot=TRUE)

dev.off()



# Save tables
normalized_dataset <- getData(methData.unite)

write.table(
      normalized_dataset,
      file = "methylkit_normalized_dataset.tsv",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

#write_methylkit_list_tables <- function(obj, prefix) {

#  dir.create(prefix, showWarnings = FALSE, recursive = TRUE)

#  for (i in seq_along(obj)) {
#    sample_name <- obj[[i]]@sample.id
##    sample_name <- gsub("[^A-Za-z0-9_.-]", "_", sample_name)

#    df <- methylKit::getData(obj[[i]])

#    write.table(
#      df,
#      file = file.path(prefix, paste0(sample_name, ".tsv")),
#      sep = "\t",
#      quote = FALSE,
#      row.names = FALSE
#    )
#  }
#}

#write_methylkit_list_tables(methData,      "tables_methData")
#write_methylkit_list_tables(methData.filt, "tables_methData_filt")
#write_methylkit_list_tables(methData.norm, "tables_methData_norm")
