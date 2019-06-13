log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tximport")

quant_files <- snakemake@input[["quant"]]
names(quant_files) <- basename(dirname(quant_files)) # dirname to drop quant.sf, basename to drop rest of path

gene2tx <- read.table(snakemake@input[["gene_trans_map"]], header=FALSE, col.names=c('gene','transcript'))
tx2gene <- gene2tx[,c('transcript', 'gene')]
## Read Salmon abundances, summarize to gene level
txi <- tximport(files = quant_files, type = "salmon", txOut = FALSE, tx2gene = tx2gene)

write.table(txi$counts, snakemake@output[["counts"]], sep = "\t", col.names = NA, quote = F)
