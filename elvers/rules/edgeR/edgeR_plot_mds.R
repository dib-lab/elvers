#from https://github.com/snakemake-workflows/rna-seq-star-deseq2
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("edgeR")

# load edgeR data
dge <- readRDS(snakemake@input[[1]])

# obtain normalized counts
#counts <- dge$counts 
# deseq2 version
#counts <- rlog(dds, blind=FALSE)


#svg(snakemake@output[[1]])
pdf(snakemake@output[[1]])
plotMDS(dge, intgroup=snakemake@params[["pca_labels"]])
dev.off()
