# Differential expression analysis with DESeq2

We will be using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to compare gene expression differences in samples between experimental conditions.

## Quickstart: Running deseq2 via eelpond

Test data:
```
./run_eelpond nema-test diffexp
```
Note that you either need to 1) have already run an assembly, such that a `fasta` file is sitting in the `eelpond/assembly` directory, 2) Run an assembly at the same time, or 3) pass an assembly in via `assemblyinput`. If you haven't run read trimming via [trimmomatic](trimmomatic.md), and read quantification via [salmon](salmon.md), `diffexp` will run these for you.

If you have not already run `./run_eelpond nema-test assemble`:

   2) Run trinity assembly at the same time:
   ```
   ./run_eelpond nema-test assemble diffexp
   ```
   3) OR, Pass an assembly in via `assemblyinput`
   ```
   ./run_eelpond assemblyinput diffexp
   ```
   with an assembly in your `yaml` configfile, e.g.:
   ```
   assemblyinput:
     assembly: rna_testdata/nema.fasta
     gene_trans_map:  rna_testdata/nema.fasta.gene_trans_map
     assembly_extension: '_input'
     ```
    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option


## DESeq2 Commands

This pipeline uses snakemake to run a few R scripts to conduct basic differential expression analysis. We read in transcript abundance information (generated with [salmon](salmon.md)) via [tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html). Note that in the salmon step, we combine files of all "units" within a sample in order to then conduct differential expression at the sample level. 

We assume the assembly has a gene-to-transcript map, such as the one produced via trinity. This is a tab separated file (`transcript \t gene`) that enables count data to be aggregated at the gene level prior to differntial expression analysis. This is recommended, see [Soneson et al, 2016](https://f1000research.com/articles/4-1521/v2). However, if you do not have this mapping, we provide an option to conduct differential expression at the transcript level via the config (see "Customizing DESeq2 Parameters" section, below).

After reading in count data, we take in two additional pieces of information: first, the sample names in the `samples.tsv` document, and second the desired `contrast`, provided as part of the DESeq2 parameters, below.

You can find these R scripts in the `eelpond` [github repo](https://github.com/dib-lab/eelpond/tree/master/rules/deseq2).

## Customizing DESeq2 parameters






## References

* [Documentation for DESeq2 with example analysis](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
* [Love et al. 2016](https://www.nature.com/nbt/journal/v34/n12/full/nbt.3682.html)

Additional links:
* [DE lecture by Jane Khudyakov, July 2017](https://github.com/Open-Data-Science-at-SIO/RNAseq-workshop-2017/blob/master/_static/Jane_differential_expression.pdf)
* [Example DE analysis from two populations of killifish! (Fundulus heteroclitus MDPL vs. MDPL)](http://htmlpreview.github.io/?https://github.com/ljcohen/Fhet_MDPL_MDPP_salinity_DE/blob/master/Fhet_MDPL_v_MDPP_interactiononly_FW_BW.html)
* [A Review of Differential Gene Expression Software for mRNA sequencing](https://github.com/ljcohen/ECE221_final_project/blob/master/Cohen_Li_ECE221_review-differential-gene.pdf)

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
