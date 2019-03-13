# Transcriptome Assembly


## Kmer Trimming 

Before running transcriptome assembly, we recommend doing some kmer spectral error trimming on your dataset, and if you have lots of reads, also performing digital normalization. We use [khmer](https://khmer.readthedocs.io/en/v2.1.1/) for both of these tasks.

You can choose whether or not to use khmer diginal normalization with `--no-diginorm`. Note that you can also conduct diginorm with the Trinity assembler.

The commands are as follows:

With digital normalizition:
```
" (interleave-reads.py {input.r1} {input.r2} && zcat {input.r1_orphan} {input.r2_orphan}) | "
" (trim-low-abund.py -V -k {params.k} -Z {params.Z} -C {params.C} - -o - -M {params.memory} "
" --diginorm --diginorm-coverage={params.cov}) | (extract-paired-reads.py --gzip "
" -p {output.paired} -s {output.single}) > {log}; split-paired-reads.py {output.paired} "
" -1 {output.r1_out} -2 {output.r2_out} >> {log}"
```

Without digital normalization:
```
" (interleave-reads.py {input.r1} {input.r2} && zcat {input.r1_orphan} {input.r2_orphan}) | "
" (trim-low-abund.py -V -k {params.k} -Z {params.Z} -C {params.C} - -o - -M {params.memory} "
"| (extract-paired-reads.py --gzip  -p {output.paired} -s {output.single}) > {log}; 
" split-paired-reads.py {output.paired} -1 {output.r1_out} -2 {output.r2_out} >> {log}"
```

## Assembling with Trinity

We use the [Trinity *de novo* transcriptome assembler](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to take short, trimmed/diginorm Illumina reads data and assemble (predict) full-length transcripts into a single fasta file output. Each contig in the fasta assembly file represents one unique transcript. The default *k*-mer size for the Trinity assembler is *k* = 25.

The resulting output assembly fasta file can then be used to align the original, trimmed (not diginorm) short Illumina reads and quantify expression per transcript.

The ID for each transcript is output (version 2.2.0 to current) as follows, where the `TRINITY` is constant, the `DN2202` is an example of a variable contig/transcript ID, `c` stands for component, `g` gene and `i` isoform:

```
TRINITY_DN2202_c0_g1_i1
```

This snakemake pipeline will run the following command:

```
Trinity --left left.fq \
  --right right.fq --seqType fq --max_memory 10G \
  --CPU 4
```

Note, the current version of Trininty (after 2.3.2) is configured to diginorm the input reads before assembly begins. Since we have already applied diginorm to our reads, the result will be a negligible decrease in read counts prior to the assembly. Applying diginorm twice is fine. For data sets with large numbers of reads, applying diginorm as a separate step as we have here may decrease the memory requirements needed by the Trinity pipeline.
