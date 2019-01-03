# Short read quality and trimming

The first step of any sequencing pipeline is to assess read quality and perform quality control as necessary. We natively enable quality assessment with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiqc](http://multiqc.info/), and do quality trimming with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

## Quality Assessment

We’re going to use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.5) and [multiqc](http://multiqc.info/) (version 1.2) to summarize the data before and after adapter trimming. There are several caveats about FastQC - the main one is that it only calculates certain statistics (like duplicated sequences) for subsets of the data (e.g. duplicate sequences are only analyzed for the first 100,000 sequences in each file.

Multiqc will summarize individual fastqc output into one output so that you can see all quality information for all files simultaenously.

## Adapter trim each pair of files

We use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.36) to trim off residual Illumina adapters that were left behind after demultiplexing.

This pipeline assumes TruSeq3-PE.fa adapters. However, if running this on your own data, you’ll need to know which Illumina sequencing adapters were used for your library prep in order to trim them off. If they are the right adapters, you should see that some of the reads are trimmed; if they’re not, you won’t see anything get trimmed.

See excellent paper on trimming parameters by [MacManes 2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full).

Based on these recommendations by MacManes 2014, we use this command in this pipeline:

```
TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:25
```


## Customizing the trimming pipeline

The default trimming paramters can be overridden by providing the following in the `.yaml` configuration file:

```
trimmomatic:
  trim_cmd:  "ILLUMINACLIP:{}:2:40:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35"
```

Be sure to modify the trimming commands as desired

