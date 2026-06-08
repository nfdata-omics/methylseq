#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(genomation)
  library(GenomicRanges)
})

option_list <- list(
  make_option("--dmr_tsv", type="character",
              help="TSV file with differentially methylated CpGs (from getData)"),
  make_option("--refseq_bed", type="character",
              help="RefSeq genes BED file"),
  make_option("--cpg_bed", type="character",
              help="CpG island BED file"),
  make_option("--out_prefix", type="character", default="annotation",
              help="Prefix for output files")
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  file.exists(opt$dmr_tsv),
  file.exists(opt$refseq_bed),
  file.exists(opt$cpg_bed)
)

dmr_tsv    <- opt$dmr_tsv
refseq_bed <- opt$refseq_bed
cpg_bed    <- opt$cpg_bed
out_prefix <- opt$out_prefix


#dmr_df     <- myDiff25p_df
#refseq_bed <- "/Users/youssef.abili/methylkit/mm10.refseq.genes.bed"
#cpg_bed    <- "/Users/youssef.abili/methylkit/mm10.cpg.bed.txt"
#out_prefix <- "out_prefix"

# -----------------------------
# Load DMCs
# -----------------------------
dmr_df <- read.delim(dmr_tsv)

#dmr_df$chr <- paste0("chr", dmr_df$chr)

dmr_gr <- makeGRangesFromDataFrame(
  dmr_df,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)

# -----------------------------
# Load gene annotation
# -----------------------------
refseq_anot <- readTranscriptFeatures(refseq_bed)

# Annotate with gene parts
dmr_gene_anot <- annotateWithGeneParts(
  target = dmr_gr,
  feature = refseq_anot
)

gene_part_df <- data.frame(
  dmr_gr,
  gene_part = getMembers(dmr_gene_anot)
)

write.table(
  gene_part_df,
  "cpg_gene_parts.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Save annotation object
#saveRDS(dmr_gene_anot, paste0(out_prefix, "_gene_annotation.rds"))

# -----------------------------
# TSS distance
# -----------------------------
dist_tss <- getAssociationWithTSS(dmr_gene_anot)

write.table(
  dist_tss,
  paste0(out_prefix, "_tss_distance.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# Gene feature summary plot
# -----------------------------
pdf(paste0(out_prefix, "_gene_parts.pdf"))

plotTargetAnnotation(dmr_gene_anot,
                     main = "Differential Methylation – Gene Annotation")

dev.off()

# -----------------------------
# CpG island / shore annotation
# -----------------------------
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

cpg_context_df <- data.frame(
  dmr_gr,
  cpg_context = getMembers(dmr_cpg_anot)
)

write.table(
  cpg_context_df,
  "cpg_context.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#saveRDS(dmr_cpg_anot, paste0(out_prefix, "_cpg_annotation.rds"))

# CpG annotation plot
pdf(paste0(out_prefix, "_cpg_context.pdf"))

plotTargetAnnotation(dmr_cpg_anot,
                     main = "Differential Methylation – CpG Context")

dev.off()

cat("Annotation completed successfully\n")