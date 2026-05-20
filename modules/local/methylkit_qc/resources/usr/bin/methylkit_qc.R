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
  make_option("--hi_perc", type="double", default=99.9),
  make_option("--destrand", type="logical", default=FALSE),
  make_option("--min_per_group", type="double", default=1)
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

#functions 
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

meta = meta[sample_ids, , drop = FALSE]
files_list = as.list(files_list)

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

meth = unite(
  methData.norm,
  destrand = destrand,
  mc.cores = cores,
  min.per.group = min_per_group
)

save(meth, file = "meth_merged_data.rda")

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

my_prcomp = PCASamples(meth, obj.return = TRUE)

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
pca_df = cbind(scores[,1:4], meta)
pca_df$sample = rownames(pca_df)

pdf("PCA_covariates.pdf")

ggplot(var_df[1:10,], aes(PC, variance)) +
  geom_col() +
  geom_text(aes(label=paste0(variance, "%")), vjust=-0.3) +
  theme_light()


for (col in colnames(meta)) {
  print(plot_pca(pca_df, "PC1", "PC2", col, var_exp))
}

dev.off()


