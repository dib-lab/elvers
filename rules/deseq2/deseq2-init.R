log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("tximport")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

quant_files <- snakemake@input[["quant"]]
names(quant_files) <- basename(dirname(quant_files)) # dirname to drop quant.sf, basename to drop rest of path

#create tx2gene from trinity gene_trans_map
gene2tx <- read.table(snakemake@input[["gene_trans_map"]], header=FALSE, col.names=c('gene','transcript'))
tx2gene <- gene2tx[,c('transcript', 'gene')]

## Read Salmon abundances
txi <- tximport(files = quant_files, type = "salmon", txOut = FALSE, tx2gene = tx2gene)

# read in sample:condition info; ensure correct ordering
sample_info <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")

print(sample_info)
print(colnames(txi$counts))

#these don't currently match, bc quant files have form: {sample}_{unit}_x_base
# either need to 1. collapse units --> samples here, or do it in salmon step. Salmon step?
stopifnot(all(sample_info$sample %in% colnames(txi$counts)))
stopifnot(all(colnames(txi$counts) %in% sample_info$sample))
sample_info <- sample_info[match(colnames(txi$counts), sample_info$sample), ]

dds <- DESeqDataSetFromMatrix(countData=txi$counts, colData=sample_info, design=~ condition)

# from rna-seq star example:
# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])
