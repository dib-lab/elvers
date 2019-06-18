log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# enable parallel processing if threads > 1
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

## Open pdf file to contain any figure generated below
#pdf(gsub("rds$", "pdf", outrds))

log <- file(snakemake@log[[1]], open="wt")

# deseq2 portion
dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])
res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
# signif results

sigRes <- subset(res, padj < 0.1)


# store results
#svg(snakemake@output[["ma_plot"]])
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]])
write.table(as.data.frame(sigRes), file=snakemake@output[["sigTable"]])
