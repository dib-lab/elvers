# Digital Normalization

In this section, we’ll apply digital normalization and variable-coverage k-mer abundance trimming to the reads prior to assembly. This has the effect of reducing the computational cost of assembly without negatively affecting the quality of the assembly.

This is all run in one command, taking the trimmed reads as input, uses the orphaned reads that survived while their mated pair did not during adapter and quality trimming. Then, low-abundance reads are trimmed to a coverage of 18 and normalized to a *k*-mer (*k* = 20) coverage of 20. 
```
(interleave-reads.py {}{}.trim_1P.fq {}{}.trim_2P.fq && zcat {}orphans.fq.gz)| \\
(trim-low-abund.py -V -k 20 -Z 18 -C 2 - -o - -M 4e9 --diginorm --diginorm-coverage=20) | \\
(extract-paired-reads.py --gzip -p {}{}.paired.gz -s {}{}.single.gz) > /dev/null
```
