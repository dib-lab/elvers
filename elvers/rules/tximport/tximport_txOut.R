log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tximport")

quant_files <- snakemake@input[["quant"]]
names(quant_files) <- basename(dirname(quant_files)) # dirname to drop quant.sf, basename to drop rest of path

## Read salmon abundances, output transcript-level count table
txi <- tximport(files = quant_files, type = "salmon", txOut = TRUE)

write.table(txi$counts, snakemake@output[["counts"]], sep = "\t", col.names = NA, quote = F)
