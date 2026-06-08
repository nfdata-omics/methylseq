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
myDiff <- calculateDiffMeth(methData.unite,
                            overdispersion = overdispersion,
                            adjust = adjust)

myDiff_df <- getData(myDiff)

# Save tables
write.table(
  myDiff_df,
  file = "myDiff_df.tsv",
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


# Barplot with positions by chromosomes (only chr with freq > 1%)
#diff_hyper = myDiff_df[myDiff_df$qvalue<0.05 & myDiff_df$meth.diff>25,]
#diff_hyper = read.delim("dma/diffMeth_hyper.tsv",h=T)
#diff_hypo = myDiff_df[myDiff_df$qvalue<0.05 & myDiff_df$meth.diff<(-25),]
#diff_hypo = read.delim("dma/diffMeth_hypo.tsv",h=T)
#diff_all = rbind(diff_hyper, diff_hypo)
#diff_all = diff_all[order(diff_all$qvalue),]
#diff_all = read.delim("dma/diffMeth_all.tsv",h=T)

chr_counts <- myDiff_df %>%
  count(chr, name="n") %>%
  mutate(freq=n/sum(n)) %>%
  filter(freq>=0.01) %>%
  arrange(desc(n))
p = ggplot(chr_counts, aes(x=reorder(chr,-n), y=n)) +
  geom_col() +
  labs(x="Chromosome", y="Count", title="Chromosome frequencies (>1% of total)") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))

chr_counts_sig <- diff_all %>%
  count(chr, name="n") %>%
  arrange(desc(n))
p1 = ggplot(chr_counts_sig, aes(x=reorder(chr,-n), y=n)) +
  geom_col() +
  labs(x="Chromosome", y="Count", title="Chromosome frequencies - Significant positions") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))

chr_counts_hyper <- diff_hyper %>%
  count(chr, name="n") %>%
  arrange(desc(n))
p2 = ggplot(chr_counts_hyper, aes(x=reorder(chr,-n), y=n)) +
  geom_col() +
  labs(x="Chromosome", y="Count", title="Chromosome frequencies - Hyper-methylated positions") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))

chr_counts_hypo <- diff_hypo %>%
  count(chr, name="n") %>%
  arrange(desc(n))
p3 = ggplot(chr_counts_hypo, aes(x=reorder(chr,-n), y=n)) +
  geom_col() +
  labs(x="Chromosome", y="Count", title="Chromosome frequencies - Hypo-methylated positions") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))

pdf("chromosome_distributions.pdf")
plot(p)
plot(p1)
plot(p2)
plot(p3)
dev.off()



cat("Differential methylation analysis completed.")
