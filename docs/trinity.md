# Trinity

The Eel Pond protocol uses the [Trinity *de novo* transcriptome assembler](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to take short, trimmed/diginorm Illumina reads data and assemble (predict) full-length transcripts into a single fasta file output. Each contig in the fasta assembly file represents one unique transcript. Trinity is a single-ksize assembler, with a default of *k* = 25.

We recommend using kmer-trimmed reads (output of khmer) as input into Triniity to reduce dataset complexity without losing valuable kmers. The resulting output assembly fasta file can then be used to align the trimmed (not diginorm) short Illumina reads and quantify expression per transcript.

Note, the current version of Trininty (after 2.3.2) is configured to diginorm the input reads before assembly begins. Since we have already applied diginorm to our reads, the result will be a negligible decrease in read counts prior to the assembly. We provide options to disable this digital normalization via the config file, but applying diginorm twice is not really a problem. For data sets with large numbers of reads, applying diginorm as a separate step as we have via khmer may decrease the memory requirements needed by the Trinity pipeline.

The ID for each transcript is output (version 2.2.0 to current) as follows, where the `TRINITY` is constant, the `DN2202` is an example of a variable contig/transcript ID, `c` stands for component, `g` gene and `i` isoform:

```
TRINITY_DN2202_c0_g1_i1
```

## Trinity Command

On the command line, the command elvers runs is approximately:
```
Trinity --left left.fq \
  --right right.fq --seqType fq --max_memory 10G \
  --CPU 4
```

But we **highly recommend** you modify `max_memory` and `CPU` to fit your data and compute resources.


## Quickstart

Run Trinity via the "default" [Eel Pond workflow](eel_pond_workflow.md) or via the [assemble subworkflow](assemble.md). To run Trinity as a standalone program, see "Advanced Usage" section below.

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Trinity will output files in the `assembly` subdirectory of this output directory. The fasta file will be `BASENAME_trinity.fasta` and the gene-trans map will be `BASENAME_trinity.fasta.gene_trans_map`. 


## Modifying Params for Trinity:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `trinity` rule, run
```
elvers config trinity --print_params
```
This will print the following:
```
  ####################  trinity  ####################
trinity:
  input_kmer_trimmed: true
  input_trimmomatic_trimmed: false
  add_single_to_paired: false # would you like to add the orphaned reads to the trinity assembly?
  max_memory: 30G
  seqtype: fq
  extra: ''
  #####################################################
```
In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra trinity parameters, e.g.:
```
  extra: '--no_normalize_reads' # to turn off Trinity's digital normalization steps 
```

Within the "default" [Eel Pond workflow](eel_pond_workflow.md) or the [assemble subworkflow](assemble.md), these options enable you to choose kmer-trimmed, quality-trimmed, or raw sequencing data as input. We recommend using kmer-trimmed reads as input. If both `input_kmer_trimmed` and `input_trimmomatic_trimmed` are `False`, we will just use raw reads from the `samples.tsv` file. 

See the [Trinity documentation](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to learn more about these parameters. Be sure the modified lines go into the config file you're using to run `elvers` (see [Understanding and Configuring Workflows](configure.md)).

## Advanced Usage: Running Trinity as a standalone rule

You can run trinity as a standalone rule, instead of withing a larger `elvers` workflow. However, to do this, you need to make sure the input files are available.

For trinity, the default input files are kmer-trimmed input data (e.g. output of khmer).

If you've already done this, you can run:
```
elvers my_config trinity
```
If not, you can run the prior steps at the same time to make sure khmer can find these input files:
```
elvers my_config get_data trimmomatic khmer trinity
```

## Snakemake rule

We wrote a [Trinity snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/trinity.html) to run Trinity.

For snakemake afficionados, see the Trinity rule on [github](https://github.com/dib-lab/elvers/blob/master/rules/trinity/trinity.rule).
