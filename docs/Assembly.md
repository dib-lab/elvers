# Assembling with Trinity

We use the [Trinity *de novo* transcriptome assembler](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.5.1) to take short, trimmed/diginorm Illumina reads data and assemble (predict) full-length transcripts into a single fasta file output. Each contig in the fasta assembly file represents one unique transcript. The default *k*-mer size for the Trinity assembler is *k* = 25.

The resulting output assembly fasta file can then be used to align the original, trimmed (not diginorm) short Illumina reads and quantify expression per transcript.

The ID for each transcript is output (version 2.2.0 to current) as follows, where the `c` stands for component, `g` gene and `i` isoform:

```
TRINITY_DN2202_c0_g1_i1
```

This snakemake pipeline will run the following command:

```
Trinity --left left.fq \
  --right right.fq --seqType fq --max_memory 10G \
  --CPU 4
```

Note, the current version of Trininty (after 2.3.2) is configured to diginorm the input reads before assembly begins. Since we have already applied diginorm to our reads, the result will be a negligible decrease in read counts prior to the assembly. Applying diginorm twice is fine. For data sets with large numbers of reads, applying diginrom as a separate step as we have here may decrease the memory requirements needed by the Trinity pipeline.
