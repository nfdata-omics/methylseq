#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(genomation)
  library(GenomicRanges)
  library(IRanges)
})

option_list <- list(
  make_option("--dmr_tsv", type = "character",
              help = "TSV file with differentially methylated CpGs/regions"),
  make_option("--caller", type = "character", default = "unknown",
              help = "Differential methylation caller: methylkit or dss"),
  make_option("--result_type", type = "character", default = "unknown",
              help = "Caller result type, such as DMC, DML, or DMR"),
  make_option("--refseq_bed", type = "character",
              help = "RefSeq genes BED file used for legacy genomation outputs"),
  make_option("--cpg_bed", type = "character",
              help = "CpG island BED file"),
  make_option("--gtf_file", type = "character", default = NULL,
              help = "GTF annotation file for transcript/gene-aware annotations"),
  make_option("--prom_up", type = "integer", default = 1000,
              help = "Promoter bases upstream of TSS [default %default]"),
  make_option("--prom_down", type = "integer", default = 200,
              help = "Promoter bases downstream of TSS [default %default]"),
  make_option("--diff_cutoff", type = "numeric", default = 25,
              help = "Methylation difference cutoff for significance labels [default %default]"),
  make_option("--qvalue_cutoff", type = "numeric", default = 0.01,
              help = "Q-value cutoff for significance labels [default %default]"),
  make_option("--out_prefix", type = "character", default = "annotation",
              help = "Prefix for output files")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  file.exists(opt$dmr_tsv),
  file.exists(opt$refseq_bed),
  file.exists(opt$cpg_bed)
)

if (is.null(opt$gtf_file) || !nzchar(opt$gtf_file)) {
  stop("--gtf_file is required")
}
if (!file.exists(opt$gtf_file)) {
  stop("gtf_file not found: ", opt$gtf_file)
}

read_bed_seqnames <- function(bed_file, nrows = 10000) {
  bed_df <- read.table(
    bed_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    nrows = nrows
  )
  unique(as.character(bed_df[[1]]))
}

harmonize_seqnames <- function(target_chr, feature_chr, feature_name) {
  target_chr <- as.character(target_chr)
  feature_chr <- unique(as.character(feature_chr))

  add_chr <- ifelse(grepl("^chr", target_chr), target_chr, paste0("chr", target_chr))
  add_chr[add_chr %in% c("chrMT", "chrM")] <- "chrM"

  remove_chr <- sub("^chr", "", target_chr)
  remove_chr[remove_chr == "M"] <- "MT"

  candidates <- list(
    original = target_chr,
    add_chr = add_chr,
    remove_chr = remove_chr
  )

  overlap_counts <- vapply(
    candidates,
    function(seqnames) length(intersect(unique(seqnames), feature_chr)),
    integer(1)
  )
  best_name <- names(which.max(overlap_counts))

  if (overlap_counts[[best_name]] == 0) {
    stop(
      "No chromosome names shared between methylation output and ", feature_name,
      ". Methylation examples: ", paste(head(unique(target_chr), 5), collapse = ", "),
      "; ", feature_name, " examples: ", paste(head(feature_chr, 5), collapse = ", ")
    )
  }

  if (best_name != "original") {
    message("Adjusted methylation chromosome names for ", feature_name, " using ", best_name)
  }

  candidates[[best_name]]
}

harmonize_granges_seqlevels <- function(gr, feature_chr, feature_name) {
  old_seqlevels <- seqlevels(gr)
  new_seqlevels <- harmonize_seqnames(old_seqlevels, feature_chr, feature_name)

  if (anyDuplicated(new_seqlevels)) {
    stop(
      "Chromosome-name harmonization produced duplicated sequence levels for ",
      feature_name
    )
  }

  seqlevels(gr) <- new_seqlevels
  gr
}

resolve_column <- function(df, candidates, required = TRUE, fallback = NULL) {
  matched <- candidates[candidates %in% colnames(df)]
  if (length(matched) > 0) {
    return(matched[[1]])
  }

  lower_names <- tolower(colnames(df))
  lower_candidates <- tolower(candidates)
  lower_match <- match(lower_candidates, lower_names, nomatch = 0)
  lower_match <- lower_match[lower_match > 0]

  if (length(lower_match) > 0) {
    return(colnames(df)[lower_match[[1]]])
  }

  if (!is.null(fallback)) {
    return(fallback)
  }

  if (required) {
    stop(
      "Cannot resolve required column. Tried: ",
      paste(candidates, collapse = ", "),
      ". Available columns: ",
      paste(colnames(df), collapse = ", ")
    )
  }

  NULL
}

collapse_hits <- function(values, query_hits, n) {
  out <- rep(NA_character_, n)
  if (length(query_hits) == 0) {
    return(out)
  }
  split_values <- split(values, query_hits)
  idx <- as.integer(names(split_values))
  out[idx] <- vapply(
    split_values,
    function(x) {
      x <- unique(na.omit(as.character(x)))
      if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
    },
    character(1)
  )
  out
}

write_table <- function(df, file) {
  write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_empty_annotation <- function(out_prefix) {
  empty_df <- data.frame()
  write_table(empty_df, paste0(out_prefix, "_cpg_gene_parts.tsv"))
  write_table(empty_df, paste0(out_prefix, "_tss_distance.tsv"))
  write_table(empty_df, paste0(out_prefix, "_cpg_context.tsv"))
  write_table(empty_df, paste0(out_prefix, "_annotated.tsv"))
  writeLines("No differential methylation loci to annotate", paste0(out_prefix, "_enrichment_signif.txt"))
  writeLines("No differential methylation loci to annotate", paste0(out_prefix, "_enrichment_direction.txt"))

  pdf(paste0(out_prefix, "_gene_parts.pdf"))
  plot.new()
  title("No differential methylation loci")
  dev.off()

  pdf(paste0(out_prefix, "_cpg_context.pdf"))
  plot.new()
  title("No differential methylation loci")
  dev.off()

  pdf(paste0(out_prefix, "_annotation_summary.pdf"))
  plot.new()
  title("No differential methylation loci")
  dev.off()

  cat("No differential methylation loci to annotate\n")
}

safe_barplot <- function(df, group_col, fill_col, main, xlab, legend_title) {
  if (!all(c(group_col, fill_col) %in% colnames(df))) {
    return(invisible(FALSE))
  }
  tab <- table(df[[group_col]], df[[fill_col]])
  if (length(tab) == 0 || any(dim(tab) == 0)) {
    return(invisible(FALSE))
  }
  barplot(
    t(tab),
    beside = FALSE,
    col = seq_len(ncol(tab)) + 1,
    main = main,
    xlab = xlab,
    ylab = "Number of loci",
    las = 2
  )
  legend("topright", legend = colnames(tab), fill = seq_len(ncol(tab)) + 1, title = legend_title, cex = 0.8)
  invisible(TRUE)
}

write_binary_enrichment <- function(df, class_col, annotation_col, label, con) {
  if (!all(c(class_col, annotation_col) %in% colnames(df))) {
    cat(label, " not available\n\n", sep = "")
    return(invisible(NULL))
  }

  tab <- table(df[[class_col]], df[[annotation_col]])
  cat("Contingency table ", annotation_col, " vs ", class_col, ":\n", sep = "")
  print(tab)
  if (all(dim(tab) == c(2, 2))) {
    cat("\nFisher's exact test:\n")
    print(fisher.test(tab))
  } else {
    cat("\nFisher's exact test skipped: expected a 2x2 table.\n")
  }
  cat("\n\n")
}

write_context_enrichment <- function(df, class_col, context_col, con) {
  if (!all(c(class_col, context_col) %in% colnames(df))) {
    cat(context_col, " not available\n\n", sep = "")
    return(invisible(NULL))
  }

  tab <- table(df[[class_col]], df[[context_col]])
  cat("Contingency table ", context_col, " vs ", class_col, ":\n", sep = "")
  print(tab)
  if (all(dim(tab) >= c(2, 2))) {
    cat("\nChi-square test:\n")
    test <- suppressWarnings(chisq.test(tab))
    print(test)
    cat("\nStandardized residuals:\n")
    print(test$stdres)
  } else {
    cat("\nChi-square test skipped: expected at least a 2x2 table.\n")
  }
  cat("\n\n")
}

read_cpg_bed <- function(cpg_bed) {
  first_line <- readLines(cpg_bed, n = 1)
  first_fields <- strsplit(first_line, "\t", fixed = TRUE)[[1]]
  has_header <- any(c("chrom", "chromStart", "chromEnd") %in% first_fields)

  if (has_header) {
    cpg_df <- read.delim(cpg_bed, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    chr_col <- resolve_column(cpg_df, c("chrom", "chr", "seqnames", "seqname"))
    start_col <- resolve_column(cpg_df, c("chromStart", "start"))
    end_col <- resolve_column(cpg_df, c("chromEnd", "end"))
    id_col <- resolve_column(cpg_df, c("ID_pos", "name", "id"), required = FALSE, fallback = NULL)
    # UCSC table/browser BED coordinates are 0-based half-open.
    starts <- as.integer(cpg_df[[start_col]]) + 1L
  } else {
    cpg_df <- read.table(cpg_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
    chr_col <- "V1"
    start_col <- "V2"
    end_col <- "V3"
    id_col <- if (ncol(cpg_df) >= 4) "V4" else NULL
    starts <- as.integer(cpg_df[[start_col]]) + 1L
  }

  ids <- if (is.null(id_col)) {
    paste0(cpg_df[[chr_col]], ":", cpg_df[[start_col]], "-", cpg_df[[end_col]])
  } else {
    as.character(cpg_df[[id_col]])
  }

  GRanges(
    seqnames = as.character(cpg_df[[chr_col]]),
    ranges = IRanges(start = starts, end = as.integer(cpg_df[[end_col]])),
    strand = "*",
    cpg_island_id = ids
  )
}

dmr_tsv <- opt$dmr_tsv
refseq_bed <- opt$refseq_bed
cpg_bed <- opt$cpg_bed
out_prefix <- opt$out_prefix
caller <- tolower(opt$caller)
result_type <- toupper(opt$result_type)

if (file.info(dmr_tsv)$size == 0) {
  write_empty_annotation(out_prefix)
  quit(save = "no", status = 0)
}

dmr_df <- read.delim(dmr_tsv, check.names = FALSE)

if (nrow(dmr_df) == 0) {
  write_empty_annotation(out_prefix)
  quit(save = "no", status = 0)
}

chr_col <- resolve_column(dmr_df, c("chr", "chrom", "chromosome", "seqnames", "seqname"))
start_col <- resolve_column(dmr_df, c("start", "pos", "position", "chr.pos", "chrPos"))
end_col <- resolve_column(dmr_df, c("end", "stop", "pos", "position", "chr.pos", "chrPos"), fallback = start_col)
strand_col <- resolve_column(dmr_df, c("strand"), required = FALSE, fallback = NULL)
pvalue_col <- resolve_column(dmr_df, c("pvalue", "p.value", "pval", "p_value"), required = FALSE, fallback = NULL)
qvalue_col <- resolve_column(dmr_df, c("qvalue", "q.value", "qval", "q_value", "fdr", "padj"), required = FALSE, fallback = NULL)
meth_diff_col <- resolve_column(
  dmr_df,
  c("meth.diff", "meth_diff", "diff.Methy", "diff", "difference", "stat"),
  required = FALSE,
  fallback = NULL
)

dmr_df$chr <- as.character(dmr_df[[chr_col]])
dmr_df$start <- as.integer(dmr_df[[start_col]])
dmr_df$end <- as.integer(dmr_df[[end_col]])
dmr_df$strand <- if (is.null(strand_col)) "*" else as.character(dmr_df[[strand_col]])

dmr_df <- dmr_df[
  !is.na(dmr_df$chr) &
    !is.na(dmr_df$start) &
    !is.na(dmr_df$end),
  ,
  drop = FALSE
]

if (nrow(dmr_df) == 0) {
  write_empty_annotation(out_prefix)
  quit(save = "no", status = 0)
}

dmr_df$start <- pmax(dmr_df$start, 1L)
swap_rows <- dmr_df$start > dmr_df$end
if (any(swap_rows)) {
  old_start <- dmr_df$start[swap_rows]
  dmr_df$start[swap_rows] <- dmr_df$end[swap_rows]
  dmr_df$end[swap_rows] <- old_start
}

refseq_chr <- read_bed_seqnames(refseq_bed)
dmr_df$chr <- harmonize_seqnames(dmr_df$chr, refseq_chr, "RefSeq BED")

dmr_gr <- makeGRangesFromDataFrame(
  dmr_df,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

annotated_df <- as.data.frame(dmr_gr)
colnames(annotated_df)[colnames(annotated_df) == "seqnames"] <- "chr"

if (!is.null(pvalue_col)) {
  annotated_df$pvalue <- dmr_df[[pvalue_col]]
}
if (!is.null(qvalue_col)) {
  annotated_df$qvalue <- dmr_df[[qvalue_col]]
}
if (!is.null(meth_diff_col)) {
  effect <- as.numeric(dmr_df[[meth_diff_col]])

  # DSS reports group1 - group2 as a fraction; this workflow defines
  # group1=control and group2=case. Normalize to case - control percentage points.
  if (caller == "dss") {
    effect <- -100 * effect
  }

  annotated_df$meth.diff <- effect
  annotated_df$meth_dir <- ifelse(annotated_df$meth.diff > 0, "hyper", "hypo")
}
if (!is.null(meth_diff_col) && !is.null(qvalue_col)) {
  annotated_df$significance <- ifelse(
    abs(annotated_df$meth.diff) >= opt$diff_cutoff & annotated_df$qvalue < opt$qvalue_cutoff,
    "signif",
    "not_signif"
  )
}

# Legacy genomation gene-part outputs retained for compatibility.
refseq_anot <- readTranscriptFeatures(refseq_bed)
dmr_gene_anot <- annotateWithGeneParts(target = dmr_gr, feature = refseq_anot)
gene_part_df <- data.frame(dmr_gr, gene_part = getMembers(dmr_gene_anot))
write_table(gene_part_df, paste0(out_prefix, "_cpg_gene_parts.tsv"))

dist_tss <- getAssociationWithTSS(dmr_gene_anot)
write_table(dist_tss, paste0(out_prefix, "_tss_distance.tsv"))

pdf(paste0(out_prefix, "_gene_parts.pdf"))
plotTargetAnnotation(dmr_gene_anot, main = "Differential Methylation - Gene Annotation")
dev.off()

# CpG island and shore annotation with strand-neutral north/south shore IDs.
cpg_islands <- read_cpg_bed(cpg_bed)
dmr_seqnames <- unique(as.character(seqnames(dmr_gr)))
cpg_islands <- cpg_islands[as.character(seqnames(cpg_islands)) %in% dmr_seqnames]
if (length(cpg_islands) == 0) {
  cpg_islands <- read_cpg_bed(cpg_bed)
  cpg_islands <- harmonize_granges_seqlevels(cpg_islands, dmr_seqnames, "methylation loci")
}

island_hits <- findOverlaps(dmr_gr, cpg_islands, ignore.strand = TRUE)
annotated_df$CpG_island_ID <- collapse_hits(mcols(cpg_islands)$cpg_island_id[subjectHits(island_hits)], queryHits(island_hits), length(dmr_gr))
annotated_df$in_CpG_island <- !is.na(annotated_df$CpG_island_ID)

north_shores <- GRanges(
  seqnames = seqnames(cpg_islands),
  ranges = IRanges(start = pmax(1L, start(cpg_islands) - 2000L), end = pmax(1L, start(cpg_islands) - 1L)),
  strand = "*",
  cpg_island_id = mcols(cpg_islands)$cpg_island_id
)
south_shores <- GRanges(
  seqnames = seqnames(cpg_islands),
  ranges = IRanges(start = end(cpg_islands) + 1L, end = end(cpg_islands) + 2000L),
  strand = "*",
  cpg_island_id = mcols(cpg_islands)$cpg_island_id
)

north_hits <- findOverlaps(dmr_gr, north_shores, ignore.strand = TRUE)
south_hits <- findOverlaps(dmr_gr, south_shores, ignore.strand = TRUE)
annotated_df$north_shore_island_ID <- collapse_hits(mcols(north_shores)$cpg_island_id[subjectHits(north_hits)], queryHits(north_hits), length(dmr_gr))
annotated_df$south_shore_island_ID <- collapse_hits(mcols(south_shores)$cpg_island_id[subjectHits(south_hits)], queryHits(south_hits), length(dmr_gr))
annotated_df$CpG_context <- "Other"
annotated_df$CpG_context[!is.na(annotated_df$north_shore_island_ID)] <- "North_shore"
annotated_df$CpG_context[!is.na(annotated_df$south_shore_island_ID)] <- "South_shore"
annotated_df$CpG_context[annotated_df$in_CpG_island] <- "CpG_island"

# Legacy CpG context output retained, but now with richer north/south IDs.
cpg_context_df <- data.frame(
  dmr_gr,
  cpg_context = annotated_df$CpG_context,
  in_CpG_island = annotated_df$in_CpG_island,
  CpG_island_ID = annotated_df$CpG_island_ID,
  north_shore_island_ID = annotated_df$north_shore_island_ID,
  south_shore_island_ID = annotated_df$south_shore_island_ID
)
write_table(cpg_context_df, paste0(out_prefix, "_cpg_context.tsv"))

# Keep the original genomation plot for broad compatibility.
cpg_anot <- readFeatureFlank(
  cpg_bed,
  feature.flank.name = c("CpGi", "shores"),
  flank = 2000
)
dmr_cpg_anot <- annotateWithFeatureFlank(
  target = dmr_gr,
  feature = cpg_anot$CpGi,
  flank = cpg_anot$shores,
  feature.name = "CpGi",
  flank.name = "shores"
)

pdf(paste0(out_prefix, "_cpg_context.pdf"))
plotTargetAnnotation(dmr_cpg_anot, main = "Differential Methylation - CpG Context")
dev.off()

# Optional GTF-driven annotations from the experimental script.
if (!is.null(opt$gtf_file)) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer is required when --gtf_file is provided")
  }

  gtf_gr <- rtracklayer::import(opt$gtf_file, format = "GTF")
  gtf_gr <- harmonize_granges_seqlevels(
    gtf_gr,
    unique(as.character(seqnames(dmr_gr))),
    "methylation loci"
  )

  transcripts_gr <- gtf_gr[gtf_gr$type == "transcript"]
  if (length(transcripts_gr) > 0) {
    keep_cols <- intersect(c("gene_id", "transcript_id", "gene_name", "gene_type", "transcript_type"), colnames(mcols(transcripts_gr)))
    transcripts_gr <- transcripts_gr[, keep_cols]

    tss_position <- ifelse(strand(transcripts_gr) == "+", start(transcripts_gr), end(transcripts_gr))
    tss_gr <- GRanges(seqnames = seqnames(transcripts_gr), ranges = IRanges(tss_position, tss_position), strand = strand(transcripts_gr))
    mcols(tss_gr) <- mcols(transcripts_gr)
    mcols(tss_gr)$TSS_position <- tss_position

    nearest_tss_hits <- distanceToNearest(dmr_gr, tss_gr, ignore.strand = TRUE, select = "all")
    nearest_tss_df <- data.frame(
      CpG_index = queryHits(nearest_tss_hits),
      TSS_index = subjectHits(nearest_tss_hits),
      stringsAsFactors = FALSE
    )
    nearest_tss_df$CpG_position <- start(dmr_gr)[nearest_tss_df$CpG_index]
    nearest_tss_df$TSS_position <- mcols(tss_gr)$TSS_position[nearest_tss_df$TSS_index]
    nearest_tss_df$TSS_strand <- as.character(strand(tss_gr))[nearest_tss_df$TSS_index]
    nearest_tss_df$distance_relative_to_TSS <- ifelse(
      nearest_tss_df$TSS_strand == "+",
      nearest_tss_df$CpG_position - nearest_tss_df$TSS_position,
      nearest_tss_df$TSS_position - nearest_tss_df$CpG_position
    )

    for (col in c("gene_name", "gene_id", "transcript_id")) {
      if (col %in% colnames(mcols(tss_gr))) {
        nearest_tss_df[[col]] <- mcols(tss_gr)[[col]][nearest_tss_df$TSS_index]
      }
    }

    promoters_gr <- promoters(transcripts_gr, upstream = opt$prom_up, downstream = opt$prom_down)
    promoter_hits <- findOverlaps(dmr_gr, promoters_gr, ignore.strand = TRUE)
    annotated_df$in_promoter <- FALSE
    annotated_df$in_promoter[unique(queryHits(promoter_hits))] <- TRUE

    if (length(promoter_hits) > 0) {
      for (col in intersect(c("gene_name", "gene_id", "transcript_id"), colnames(mcols(promoters_gr)))) {
        annotated_df[[paste0("promoter_", col)]] <- collapse_hits(
          mcols(promoters_gr)[[col]][subjectHits(promoter_hits)],
          queryHits(promoter_hits),
          length(dmr_gr)
        )
      }
    }

    annotated_df$nearest_TSS_distance <- collapse_hits(nearest_tss_df$distance_relative_to_TSS, nearest_tss_df$CpG_index, length(dmr_gr))
    for (col in intersect(c("gene_name", "gene_id", "transcript_id"), colnames(nearest_tss_df))) {
      annotated_df[[paste0("nearest_TSS_", col)]] <- collapse_hits(nearest_tss_df[[col]], nearest_tss_df$CpG_index, length(dmr_gr))
    }
  }

  genes_gr <- gtf_gr[gtf_gr$type == "gene"]
  if (length(genes_gr) > 0) {
    gene_hits <- findOverlaps(dmr_gr, genes_gr, ignore.strand = TRUE)
    annotated_df$in_gene_body <- FALSE
    annotated_df$in_gene_body[unique(queryHits(gene_hits))] <- TRUE

    for (col in intersect(c("gene_name", "gene_id", "gene_type"), colnames(mcols(genes_gr)))) {
      annotated_df[[paste0("gene_body_", col)]] <- collapse_hits(
        mcols(genes_gr)[[col]][subjectHits(gene_hits)],
        queryHits(gene_hits),
        length(dmr_gr)
      )
    }
  }
} else {
  annotated_df$gene_part <- gene_part_df$gene_part
}

write_table(annotated_df, paste0(out_prefix, "_annotated.tsv"))

pdf(paste0(out_prefix, "_annotation_summary.pdf"))
plots_written <- 0L
if ("meth_dir" %in% colnames(annotated_df)) {
  plots_written <- plots_written + as.integer(safe_barplot(annotated_df, "CpG_context", "meth_dir", "CpG Context by Methylation Direction", "CpG context", "Direction"))
}
if ("significance" %in% colnames(annotated_df)) {
  plots_written <- plots_written + as.integer(safe_barplot(annotated_df, "CpG_context", "significance", "CpG Context by Significance", "CpG context", "Significance"))
}
if (all(c("in_promoter", "meth_dir") %in% colnames(annotated_df))) {
  plots_written <- plots_written + as.integer(safe_barplot(annotated_df, "in_promoter", "meth_dir", "Promoter Overlap by Methylation Direction", "In promoter", "Direction"))
}
if (all(c("in_gene_body", "meth_dir") %in% colnames(annotated_df))) {
  plots_written <- plots_written + as.integer(safe_barplot(annotated_df, "in_gene_body", "meth_dir", "Gene Body Overlap by Methylation Direction", "In gene body", "Direction"))
}
if (plots_written == 0L) {
  plot.new()
  title("No annotation summary plots available")
  text(0.5, 0.45, paste("Caller:", caller, "Result type:", result_type))
}
dev.off()

capture.output({
  cat("SIGNIFICANCE ENRICHMENT\n\n")
  write_binary_enrichment(annotated_df, "significance", "in_promoter", "Promoter annotation", stdout())
  write_binary_enrichment(annotated_df, "significance", "in_gene_body", "Gene body annotation", stdout())
  write_context_enrichment(annotated_df, "significance", "CpG_context", stdout())
}, file = paste0(out_prefix, "_enrichment_signif.txt"))

capture.output({
  cat("METHYLATION DIRECTION ENRICHMENT\n\n")
  cat("Direction is normalized as case minus control: hyper > 0, hypo <= 0.\n")
  if ("meth_dir" %in% colnames(annotated_df)) {
    cat("Caller: ", caller, "; result type: ", result_type, "\n\n", sep = "")
    cat("Direction counts:\n")
    print(table(annotated_df$meth_dir))
    cat("\n")
  }
  write_binary_enrichment(annotated_df, "meth_dir", "in_promoter", "Promoter annotation", stdout())
  write_binary_enrichment(annotated_df, "meth_dir", "in_gene_body", "Gene body annotation", stdout())
  write_context_enrichment(annotated_df, "meth_dir", "CpG_context", stdout())
}, file = paste0(out_prefix, "_enrichment_direction.txt"))

cat("Annotation completed successfully\n")
