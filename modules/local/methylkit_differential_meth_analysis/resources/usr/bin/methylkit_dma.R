#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(methylKit)
})

option_list = list(
  make_option("--meth_rda", type="character", help="Path to merged meth object RDA file from previous step"),
  #make_option("--output_dir", type="character", help="Directory to save differential methylation results"),
  make_option("--diff_cutoff", type="numeric", default=25, help="Methylation difference cutoff"),
  make_option("--qvalue_cutoff", type="numeric", default=0.01, help="Q-value cutoff for significance"),
  make_option("--overdispersion", type="character", default="MN", help="Overdispersion model for calculateDiffMeth"),
  make_option("--adjust", type="character", default="BH", help="Multiple testing correction method")
)

opt = parse_args(OptionParser(option_list=option_list))

meth_rda <- opt$meth_rda
#output_dir <- opt$output_dir
diff_cutoff <- opt$diff_cutoff
qvalue_cutoff <- opt$qvalue_cutoff
overdispersion <- opt$overdispersion
adjust <- opt$adjust

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
load(meth_rda) # assumes 'meth' object is loaded

# Fix absolute DB paths
dbpath <- meth@dbpath 
dbpath <- paste( "./methylDB_dir/methylDB/" , basename(dbpath), sep="")
meth@dbpath <- dbpath

# Differential methylation
myDiff <- calculateDiffMeth(meth,
                            overdispersion = overdispersion,
                            adjust = adjust)

myDiff_df <- getData(myDiff)

# Volcano plot
pdf( "volcano_plot.pdf")
plot(myDiff_df$meth.diff, -log10(myDiff_df$qvalue),
     xlab = "Methylation difference (%)",
     ylab = "-log10(Q-value)",
     main = "Differential Methylation Volcano Plot")
abline(v=0, col="red")
dev.off()

# Differential methylation per chromosome
pdf( "diffMeth_per_chr.pdf")
diffMethPerChr(myDiff)
dev.off()

# Extract hyper, hypo, and all significant DMCs
myDiff.hyper <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = qvalue_cutoff, type = "hyper")
myDiff.hypo  <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = qvalue_cutoff, type = "hypo")
myDiff.all   <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = qvalue_cutoff)

# Order by qvalue
myDiff25p.hyper_df <- getData(myDiff.hyper)
myDiff25p.hyper_df <- myDiff25p.hyper_df[order(myDiff25p.hyper_df$qvalue), ]

myDiff25p.hypo_df <- getData(myDiff.hypo)
myDiff25p.hypo_df  <- myDiff25p.hypo_df[order(myDiff25p.hypo_df$qvalue), ]

myDiff.all_df <- getData(myDiff.all)
myDiff.all_df   <- myDiff.all_df[order(myDiff.all_df$qvalue), ]

# Save results as TSV
write.table(myDiff25p.hyper_df, file =  "diffMeth_hyper.tsv", sep="\t", row.names = FALSE, quote = FALSE)
write.table(myDiff25p.hypo_df,  file =  "diffMeth_hypo.tsv", sep="\t", row.names = FALSE, quote = FALSE)
write.table(myDiff.all_df,   file =  "diffMeth_all.tsv", sep="\t", row.names = FALSE, quote = FALSE)

cat("Differential methylation analysis completed.")
