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

## Read Salmon abundances
# when create tx2gene, make sure cols are in appropriate order & have appropriate headers
txi <- tximport(files = quant_files, type = "salmon", txOut = FALSE, tx2gene = snakemake@input[["tx2gene"]])

# read in sample:condition info; ensure correct ordering
sample_info <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")
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
