log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("tximport")
library("stringr")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

quant_files <- snakemake@input[["quant"]]
names(quant_files) <- basename(dirname(quant_files)) # dirname to drop quant.sf, basename to drop rest of path

init_gtm <- function(gtm, quant_files){
    gene2tx <- read.table(gtm, header=FALSE, col.names=c('gene','transcript'))
    tx2gene <- gene2tx[,c('transcript', 'gene')]
    ## Read Salmon abundances, summarize to gene level
    txi <- tximport(files = quant_files, type = "salmon", txOut = FALSE, tx2gene = tx2gene)
    return(txi)
}

init_txout <- function(qfiles){
    ## Read Salmon abundances, transcript level
    txi <- tximport(files = qfiles, type = "salmon", txOut = TRUE)
    print("txOut")
    print(txi)
    return(txi)
}

txi <- tryCatch(
    {
    init_gtm(snakemake@input[["gene_trans_map"]], quant_files)
    },
    error = function(e){
        init_txout(quant_files)
    }
)

# read in sample:condition info; ensure correct ordering
sample_info <- read.delim(snakemake@params[["samples"]], header=TRUE) 

# remove excess from sample names
colnames(txi$counts) <- str_extract(colnames(txi$counts), "[^_]+")

# check all samples are present in both
stopifnot(all(sample_info$sample %in% colnames(txi$counts)))
stopifnot(all(colnames(txi$counts) %in% sample_info$sample))
# reorder sample info to match tximport
sample_info <- sample_info[match(colnames(txi$counts),sample_info$sample), ]
row.names(sample_info) <- sample_info$sample
sample_info$sample <- NULL
sample_info$condition <- factor(sample_info$condition)

# generate DESeq data set (dds)
dds <- DESeqDataSetFromTximport(txi, sample_info, ~condition)

#pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10 # just reduce dataset size. More stringent independent filtering happens later within DESeq function.
dds <- dds[keep,]

# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[[1]])
