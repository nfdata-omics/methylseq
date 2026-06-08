#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DSS)
  library(bsseq)
})

option_list <- list(
  make_option("--cov_dir", type = "character", help = "Directory with Bismark .cov files"),
  make_option("--metadata", type = "character", help = "Metadata CSV/TSV file"),
  make_option("--sample_suffix", type = "character", default = "_bismark.cov.gz"),
  make_option("--pattern", type = "character", default = "\\.cov"),
  make_option("--group_column", type = "character", help = "Metadata column defining groups"),
  make_option("--group1", type = "character", help = "First group / condition"),
  make_option("--group2", type = "character", help = "Second group / condition"),
  make_option("--sep", type = "character", default = ",", help = "Metadata separator"),
  make_option("--min_coverage", type = "integer", default = 5),
  make_option("--smoothing_single", type = "logical", default = TRUE),
  make_option("--smoothing_replicates", type = "logical", default = FALSE),
  make_option("--smoothing_span", type = "integer", default = NA),
  make_option("--p_threshold", type = "double", default = 0.001),
  make_option("--delta", type = "double", default = 0.10),
  make_option("--minlen", type = "integer", default = 20),
  make_option("--minCG", type = "integer", default = 3),
  make_option("--dis_merge", type = "integer", default = 50),
  make_option("--out_prefix", type = "character", default = NA)
)

opt <- parse_args(OptionParser(option_list = option_list))

cov_dir <- opt$cov_dir
metadata <- opt$metadata
sample_suffix <- opt$sample_suffix
pattern <- opt$pattern
group_column <- opt$group_column
group1 <- opt$group1
group2 <- opt$group2
sep <- opt$sep

#debug
#cov_dir <- "/Users/youssef.abili/HT/DATA/sallese/dss-single/input"
#metadata <- "/Users/youssef.abili/HT/DATA/sallese/dss-single/input/metadata.csv"
#sample_suffix <- "_1_val_1_bismark_bt2_pe.bismark.cov.gz"
#pattern <- "\\.bismark\\.cov\\.gz$"
#group_column <- "genotype"
#group1 <- "WT"
#group2 <- "Mut"
#sep <- ","

out_prefix <- if (is.na(opt$out_prefix)) {
  paste0(group1, "_vs_", group2)
} else {
  opt$out_prefix
}

smoothing_span <- if (is.na(opt$smoothing_span)) NULL else opt$smoothing_span

stopifnot(
  dir.exists(cov_dir),
  file.exists(metadata),
  !is.null(group_column),
  !is.null(group1),
  !is.null(group2)
)

files <- list.files(
  path = cov_dir,
  pattern = pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No coverage files found in: ", cov_dir)
}

sample_ids <- sub(sample_suffix, "", basename(files), fixed = TRUE)

file_table <- data.frame(
  sample = sample_ids,
  file = files,
  stringsAsFactors = FALSE
)

meta <- read.delim(
  metadata,
  row.names = 1,
  check.names = FALSE,
  sep = sep
)

if (!group_column %in% colnames(meta)) {
  stop("Group column not found in metadata: ", group_column)
}

meta <- meta[rownames(meta) %in% file_table$sample, , drop = FALSE]
file_table <- file_table[file_table$sample %in% rownames(meta), , drop = FALSE]

if (nrow(meta) == 0) {
  stop("No matching samples between metadata and coverage files.")
}

group1_samples <- rownames(meta)[meta[[group_column]] == group1]
group2_samples <- rownames(meta)[meta[[group_column]] == group2]

if (length(group1_samples) == 0) {
  stop("No samples found for group1: ", group1)
}

if (length(group2_samples) == 0) {
  stop("No samples found for group2: ", group2)
}

selected_samples <- c(group1_samples, group2_samples)

file_table <- file_table[
  match(selected_samples, file_table$sample),
  ,
  drop = FALSE
]

if (any(is.na(file_table$file))) {
  stop("Missing coverage files for selected samples.")
}

n_group1 <- length(group1_samples)
n_group2 <- length(group2_samples)

if (n_group1 == 1 && n_group2 == 1) {
  dss_mode <- "DSS-single"
  smoothing <- opt$smoothing_single
} else if (n_group1 >= 2 && n_group2 >= 2) {
  dss_mode <- "DSS-replicate"
  smoothing <- opt$smoothing_replicates
} else {
  stop(
    "Unsupported DSS design: ",
    group1, " has ", n_group1, " sample(s), ",
    group2, " has ", n_group2, " sample(s). ",
    "Use either 1 vs 1 or replicated groups on both sides."
  )
}

message("Running mode: ", dss_mode)
message(group1, ": ", paste(group1_samples, collapse = ", "))
message(group2, ": ", paste(group2_samples, collapse = ", "))

convert_bismark_cov <- function(file, min_coverage = 5) {
  dt <- fread(file)
  
  if (ncol(dt) < 6) {
    stop("Unexpected Bismark coverage format in file: ", file)
  }
  
  dt <- dt[, 1:6]
  setnames(dt, c("chr", "start", "end", "pct", "meth", "unmeth"))
  
  dt[, chr := as.character(chr)]
  dt[, pos := as.integer(start)]
  dt[, X := as.integer(meth)]
  dt[, N := as.integer(meth + unmeth)]
  
  dt <- dt[!is.na(chr) & !is.na(pos) & !is.na(N) & !is.na(X)]
  dt <- dt[N >= min_coverage]
  dt <- dt[X <= N]
  
  dt[, .(chr, pos, N, X)]
}

bs_list <- lapply(
  file_table$file,
  convert_bismark_cov,
  min_coverage = opt$min_coverage
)

names(bs_list) <- file_table$sample

bsobj <- makeBSseqData(
  dat = bs_list,
  sampleNames = file_table$sample
)

dml_args <- list(
  BSobj = bsobj,
  group1 = group1_samples,
  group2 = group2_samples,
  smoothing = smoothing
)

if (!is.null(smoothing_span)) {
  dml_args$smoothing.span <- smoothing_span
}

dml_test <- do.call(DMLtest, dml_args)

dmrs <- callDMR(
  dml_test,
  p.threshold = opt$p_threshold,
  delta = opt$delta,
  minlen = opt$minlen,
  minCG = opt$minCG,
  dis.merge = opt$dis_merge
)

dml_file <- paste0(out_prefix, "_", dss_mode, "_DML.tsv")
dmr_file <- paste0(out_prefix, "_", dss_mode, "_DMR.tsv")
rda_file <- paste0(out_prefix, "_", dss_mode, "_objects.rda")
summary_file <- paste0(out_prefix, "_", dss_mode, "_summary.tsv")

fwrite(
  as.data.table(dml_test),
  dml_file,
  sep = "\t"
)

if (!is.null(dmrs) && nrow(dmrs) > 0) {
  fwrite(
    as.data.table(dmrs),
    dmr_file,
    sep = "\t"
  )
} else {
  fwrite(
    data.table(),
    dmr_file,
    sep = "\t"
  )
  message("No DMRs found for ", group1, " vs ", group2)
}

run_summary <- data.table(
  comparison = paste0(group1, "_vs_", group2),
  dss_mode = dss_mode,
  group1 = group1,
  group2 = group2,
  n_group1 = n_group1,
  n_group2 = n_group2,
  group1_samples = paste(group1_samples, collapse = ","),
  group2_samples = paste(group2_samples, collapse = ","),
  smoothing = smoothing,
  min_coverage = opt$min_coverage,
  p_threshold = opt$p_threshold,
  delta = opt$delta,
  minlen = opt$minlen,
  minCG = opt$minCG,
  dis_merge = opt$dis_merge
)

fwrite(run_summary, summary_file, sep = "\t")

save(
  bsobj,
  dml_test,
  dmrs,
  meta,
  file_table,
  run_summary,
  file = rda_file
)