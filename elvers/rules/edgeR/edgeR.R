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
# TODO explore IHW usage


# store results
#svg(snakemake@output[["ma_plot"]])
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]])


#init: https://github.com/csoneson/rnaseqworkflow/blob/master/scripts/run_dge_edgeR.R

## Define design. ************** MODIFY ************** 
stopifnot(all(colnames(dge0) == metadata$ID))
(des <- model.matrix(~ XXXX, data = metadata))

## Filter out genes with average CPM below 1
print(dim(dge0))
cpms <- cpm(dge0)
dge <- dge0[apply(cpms, 1, mean) > 1, ]
dge <- calcNormFactors(dge)
print(dim(dge))

## Add gene annotation
annot <- tx2gene %>% dplyr::select(-tx, -tx_biotype, -start, -end) %>% distinct()
if (any(duplicated(annot$gene))) {
  stop(paste0("The following genes are represented by multiple rows in the ", 
              "gene annotation: ", 
              paste(annot$gene[duplicated(annot$gene)], collapse = ",")))
}
annot <- annot[match(rownames(dge), annot$gene), ]
rownames(annot) <- annot$gene
dge$genes <- annot

## Estimate dispersion and fit model
dge <- estimateDisp(dge, design = des)
qlfit <- glmQLFit(dge, design = des)

## Plot dispersions
plotBCV(dge)

## Define contrasts. ************** MODIFY ************** 
(contrasts <- as.data.frame(makeContrasts(XXXX, levels = des)))

## Perform tests
signif3 <- function(x) signif(x, digits = 3)
edgeR_res <- lapply(contrasts, function(cm) {
  qlf <- glmQLFTest(qlfit, contrast = cm)
  tt <- topTags(qlf, n = Inf, sort.by = "none")$table
  tt %>% dplyr::mutate_if(is.numeric, signif3)
})

## Write results to text files and make MA plots
if (class(edgeR_res) == "data.frame") {
  write.table(edgeR_res %>% dplyr::arrange(PValue), 
              file = gsub("rds$", "txt", outrds), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  print(ggplot(edgeR_res, aes(x = logCPM, y = logFC, color = FDR <= 0.05)) + 
          geom_point() + theme_bw() + 
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")))
} else {
  for (nm in names(edgeR_res)) {
    write.table(edgeR_res[[nm]] %>% dplyr::arrange(PValue), 
                file = gsub("\\.rds$", paste0("_", nm, ".txt"), outrds), 
                sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    print(ggplot(edgeR_res[[nm]], aes(x = logCPM, y = logFC, color = FDR <= 0.05)) + 
            geom_point() + theme_bw() + 
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + 
            ggtitle(nm))
  }
}

dev.off()

saveRDS(list(results = edgeR_res, data = dge0), file = outrds)

sessionInfo()
date()


