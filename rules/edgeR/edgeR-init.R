log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("edgeR")
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

#create tx2gene from trinity gene_trans_map
gene2tx <- read.table(snakemake@input[["gene_trans_map"]], header=FALSE, col.names=c('gene','transcript'))
tx2gene <- gene2tx[,c('transcript', 'gene')]

## Read Salmon abundances
txi <- tximport(files = quant_files, type = "salmon", txOut = FALSE, tx2gene = tx2gene)

# read in sample:condition info; ensure correct ordering
sample_info <- read.table(snakemake@params[["samples"]], header=TRUE)

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

##### edgeR-specific portion #####

# generate edgeR data set 
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#edger
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide

# edgeR dataset (dge)
dge <- DGEList(y, group =sample_info$condition)  #does parallel work here?
# normalization and preprocessing
dge <- calcNormFactors(dge)

# this isn't working properly -- might just be the lack of replicates, though...
design <- model.matrix(~sample_info$condition)

#dge <- estimateDisp(dge, sample_info$condition)
dge <- estimateDisp(dge)

saveRDS(dge, file=snakemake@output[[1]])
