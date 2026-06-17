#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(methylKit)
  library(ggplot2)
})

option_list = list(
  make_option("--meth_rda", type="character", help="Path to merged meth object RDA file from previous step"),
  #make_option("--output_dir", type="character", help="Directory to save differential methylation results"),
  make_option("--diff_cutoff", type="numeric", default=25, help="Methylation difference cutoff"),
  make_option("--qvalue_cutoff", type="numeric", default=0.01, help="Q-value cutoff for significance"),
  make_option("--overdispersion", type="character", default="MN", help="Overdispersion model for calculateDiffMeth"),
  make_option("--adjust", type="character", default="BH", help="Multiple testing correction method"),
  make_option("--test", type="character", default="Chisq", help="Statistical test for calculateDiffMeth"),
  make_option("--cores", type="integer", default=1, help="Number of cores for methylKit"),
  make_option("--comparison_id", type="character"),
  make_option("--case_samples", type="character"),
  make_option("--control_samples", type="character")
)

opt = parse_args(OptionParser(option_list=option_list))

meth_rda <- opt$meth_rda
#output_dir <- opt$output_dir
diff_cutoff <- opt$diff_cutoff
qvalue_cutoff <- opt$qvalue_cutoff
overdispersion <- opt$overdispersion
adjust <- opt$adjust
test <- opt$test
cores <- opt$cores

overdispersion <- match.arg(overdispersion, c("none", "MN", "shrinkMN"))
adjust <- match.arg(adjust, c("SLIM", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none", "qvalue"))
test <- match.arg(test, c("F", "Chisq", "fast.fisher", "midPval"))

if (test == "F" && overdispersion %in% c("MN", "shrinkMN")) {
  message(
    "Using test='F' with overdispersion='", overdispersion,
    "' can fail in methylKit when a locus has too few non-missing samples. ",
    "Use --test Chisq if this run hits 'argument is of length zero'."
  )
}

#DEBUG
#meth_rda <- opt$meth_rda
#output_dir <- opt$output_dir
#diff_cutoff <- 25
#qvalue_cutoff <- 0.01
#overdispersion <- "NM"
#adjust <- "BH"

stopifnot(file.exists(meth_rda))
#if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# Load merged methylation data
load(meth_rda) # assumes 'methData.unite' object is loaded

if (!exists("methData.unite")) {
  stop("Expected object 'methData.unite' was not found in ", meth_rda)
}

# Fix absolute DB paths
dbpath <- file.path("methylDB_dir", "methylDB", basename(methData.unite@dbpath))
if (!file.exists(dbpath)) {
  stop("Expected staged methylKit DB file was not found: ", dbpath)
}
methData.unite@dbpath <- dbpath


############################DEV-START

comparison_id <- opt$comparison_id
case_samples <- strsplit(opt$case_samples, ",")[[1]]
control_samples <- strsplit(opt$control_samples, ",")[[1]]

wanted_samples <- c(control_samples, case_samples)

if (!exists("meta")) {
  stop("Expected object 'meta' was not found in ", meth_rda)
}

if (!all(wanted_samples %in% rownames(meta))) {
  stop(
    "Some comparison samples are missing from metadata: ",
    paste(setdiff(wanted_samples, rownames(meta)), collapse = ", ")
  )
}

sample_order <- rownames(meta)
comparison_sample_ids <- sample_order[sample_order %in% wanted_samples]

comparison_treatment <- ifelse(
  comparison_sample_ids %in% case_samples,
  1,
  0
)

methData.unite <- reorganize(
  methData.unite,
  sample.ids = comparison_sample_ids,
  treatment = comparison_treatment
)

message("Running comparison: ", comparison_id)
message("Case samples: ", paste(case_samples, collapse = ", "))
message("Control samples: ", paste(control_samples, collapse = ", "))

############################DEV-END


# Differential methylation
myDiff <- calculateDiffMeth(methData.unite,
                            overdispersion = overdispersion,
                            adjust = adjust,
                            test = test,
                            mc.cores = cores)

myDiff_df <- getData(myDiff)

# Save tables
write.table(
  myDiff_df,
  file = paste0(comparison_id, "_diffmeth_raw.tsv") ,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Volcano plot
#pdf( "volcano_plot.pdf")
#plot(myDiff_df$meth.diff, -log10(myDiff_df$qvalue),
#     xlab = "Methylation difference (%)",
#     ylab = "-log10(Q-value)",
#     main = "Differential Methylation Volcano Plot")
#abline(v=0, col="red")
#dev.off()

# Differential methylation per chromosome
pdf( paste0(comparison_id, "_diffmeth_per_chr.pdf") )
diffMethPerChr(myDiff)
dev.off()

# Extract hyper, hypo, and all significant DMCs
myDiff.hyper <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = qvalue_cutoff, type = "hyper")
myDiff.hypo  <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = qvalue_cutoff, type = "hypo")
myDiff.all   <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = qvalue_cutoff)

# Order by qvalue
order_by_qvalue <- function(df) {
  if (!"qvalue" %in% colnames(df) || nrow(df) == 0) {
    return(df)
  }
  df[order(df$qvalue), , drop = FALSE]
}

myDiff25p.hyper_df <- order_by_qvalue(getData(myDiff.hyper))
myDiff25p.hypo_df <- order_by_qvalue(getData(myDiff.hypo))
myDiff.all_df <- order_by_qvalue(getData(myDiff.all))

# Save results as TSV
write.table(myDiff25p.hyper_df, file =  paste0(comparison_id,"_diffmeth_hyper.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
write.table(myDiff25p.hypo_df,  file =  paste0(comparison_id,"_diffmeth_hypo.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
write.table(myDiff.all_df,   file =  paste0(comparison_id,"_diffmeth_all.tsv"), sep="\t", row.names = FALSE, quote = FALSE)


chromosome_counts <- function(df, min_freq = NULL) {
  if (nrow(df) == 0 || !"chr" %in% colnames(df)) {
    return(data.frame(chr = character(), n = integer(), freq = numeric()))
  }

  counts <- as.data.frame(table(df$chr), stringsAsFactors = FALSE)
  colnames(counts) <- c("chr", "n")
  counts <- counts[counts$n > 0, , drop = FALSE]
  counts$freq <- counts$n / sum(counts$n)

  if (!is.null(min_freq)) {
    counts <- counts[counts$freq >= min_freq, , drop = FALSE]
  }

  counts[order(-counts$n), , drop = FALSE]
}

plot_chr_counts <- function(counts, title) {
  if (nrow(counts) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No loci") +
        labs(x = "Chromosome", y = "Count", title = title) +
        theme_minimal() +
        theme(
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()
        )
    )
  }

  ggplot(counts, aes(x = reorder(chr, -n), y = n)) +
    geom_col() +
    labs(x = "Chromosome", y = "Count", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Barplot with positions by chromosomes (only chr with freq > 1%)
chr_counts <- chromosome_counts(myDiff_df, min_freq = 0.01)
p <- plot_chr_counts(chr_counts, "Chromosome frequencies (>1% of total)")

chr_counts_sig <- chromosome_counts(myDiff.all_df)
p1 <- plot_chr_counts(chr_counts_sig, "Chromosome frequencies - Significant positions")

chr_counts_hyper <- chromosome_counts(myDiff25p.hyper_df)
p2 <- plot_chr_counts(chr_counts_hyper, "Chromosome frequencies - Hyper-methylated positions")

chr_counts_hypo <- chromosome_counts(myDiff25p.hypo_df)
p3 <- plot_chr_counts(chr_counts_hypo, "Chromosome frequencies - Hypo-methylated positions")

pdf(paste0(comparison_id, "_chromosome_distributions.pdf"))
plot(p)
plot(p1)
plot(p2)
plot(p3)
dev.off()



cat("Differential methylation analysis completed.")
