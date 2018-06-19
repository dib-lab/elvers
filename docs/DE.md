# Differential expression analysis with DESeq2

Comparing gene expression differences in samples between experimental conditions. 

We will be using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

References:
* [Documentation for DESeq2 with example analysis](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
* [Love et al. 2016](https://www.nature.com/nbt/journal/v34/n12/full/nbt.3682.html)

Additional links:
* [DE lecture by Jane Khudyakov, July 2017](_static/Jane_differential_expression.pdf)
* [Example DE analysis from two populations of killifish! (Fundulus heteroclitus MDPL vs. MDPL)](http://htmlpreview.github.io/?https://github.com/ljcohen/Fhet_MDPL_MDPP_salinity_DE/blob/master/Fhet_MDPL_v_MDPP_interactiononly_FW_BW.html)

## RStudio!

The pipeline will be running these commands in an R script. You could run them in R Studio:

Load libraries
```
library(DESeq2)
library("lattice")
library(tximport)
library(readr)
library(gplots)
library(RColorBrewer)
source('~/plotPCAWithSampleNames.R')
```

Tell RStudio where your files are and ask whether they exist:
```
setwd("/mnt/work/quant/salmon_out/")
dir<-"/mnt/work/quant/"
files_list = list.files()
files <- file.path(dir, "salmon_out",files_list, "quant.sf")
names(files) <- c("0Hour_1","0Hour_2","0Hour_3","0Hour_4","0Hour_5","6Hour_1","6Hour_2","6Hour_3","6Hour_4","6Hour_5")
files
print(file.exists(files))
```

Grab the [gene names](https://raw.githubusercontent.com/Open-Data-Science-at-SIO/RNAseq-workshop-2017/master/_static/nema_transcript_gene_id.txt) and transcript ID file to [summarize expression at the gene level](https://f1000research.com/articles/4-1521/v2).

```
tx2gene <- read.table("~/nema_transcript_gene_id.txt",sep="\t")
cols<-c("transcript_id","gene_id")
colnames(tx2gene)<-cols
head(tx2gene)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,importer=read.delim)
head(txi.salmon$counts)
dim(txi.salmon$counts)
```
Assign experimental variables:

```
condition = factor(c("0Hour","0Hour","0Hour","0Hour","0Hour","6Hour","6Hour","6Hour","6Hour","6Hour"))
ExpDesign <- data.frame(row.names=colnames(txi.salmon$counts), condition = condition)
ExpDesign
```

Run DESeq2:

```
dds <- DESeqDataSetFromTximport(txi.salmon, ExpDesign, ~condition)
dds <- DESeq(dds, betaPrior=FALSE)
```

Get counts:
```
counts_table = counts( dds, normalized=TRUE )
```

Filtering out low expression transcripts:

See plot from [Lisa Komoroske](https://github.com/Open-Data-Science-at-SIO/RNAseq-workshop-2017/blob/master/_static/Before-after_filter.pdf) generated with [RNAseq123](https://www.bioconductor.org/help/workflows/RNAseq123/)
```
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
head(filtered_norm_counts)
```

Estimate dispersion:

```
plotDispEsts(dds)
```

PCA:
```
log_dds<-rlog(dds)
plotPCAWithSampleNames(log_dds, intgroup="condition", ntop=40000)
```

Get DE results:

```
res<-results(dds,contrast=c("condition","6Hour","0Hour"))
head(res)
res_ordered<-res[order(res$padj),]
GeneID<-rownames(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_genes<-cbind(res_ordered,GeneID)
dim(res_genes)
head(res_genes)
dim(res_genes)
res_genes_merged <- merge(res_genes,filtered_norm_counts,by=unique("GeneID"))
dim(res_genes_merged)
head(res_genes_merged)
res_ordered<-res_genes_merged[order(res_genes_merged$padj),]
write.csv(res_ordered, file="nema_DESeq_all.csv" )
```

Set a threshold cutoff of padj<0.05 and ± log2FC 1:

```
resSig = res_ordered[res_ordered$padj < 0.05, ]
resSig = resSig[resSig$log2FoldChange > 1 | resSig$log2FoldChange < -1,]
write.csv(resSig,file="nema_DESeq_padj0.05_log2FC1.csv")
```


MA plot with gene names:

```
plot(log2(res_ordered$baseMean), res_ordered$log2FoldChange, col=ifelse(res_ordered$padj < 0.05, "red","gray67"),main="nema (padj<0.05, log2FC = ±1)",xlim=c(1,20),pch=20,cex=1,ylim=c(-12,12))
abline(h=c(-1,1), col="blue")
genes<-resSig$GeneID
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=genes,pos=2,cex=0.60)
```

Heatmap

```
d<-resSig
dim(d)
head(d)
colnames(d)
d<-d[,c(8:17)]
d<-as.matrix(d)
d<-as.data.frame(d)
d<-as.matrix(d)
rownames(d) <- resSig[,1]
head(d)

hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="nema (padj<0.05, log2FC = ±1)", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
```
