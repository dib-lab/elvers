# Trinity

## Quickstart: running trinity via eelpond:

```
./run_eelpond nema-test assemble
```
This will run preprocessing and kmer-trimming for you prior to assembly. For more options, read below!

## Trinity

The Eel Pond protocol uses  the [Trinity *de novo* transcriptome assembler](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to take short, trimmed/diginorm Illumina reads data and assemble (predict) full-length transcripts into a single fasta file output. Each contig in the fasta assembly file represents one unique transcript. Trinity is a single-ksize assembler, with a default of *k* = 25.

We recommend using kmer-trimmed reads (output of khmer) as input into Triniity to reduce dataset complexity without losing valuable kmers. The resulting output assembly fasta file can then be used to align the trimmed (not diginorm) short Illumina reads and quantify expression per transcript.

Note, the current version of Trininty (after 2.3.2) is configured to diginorm the input reads before assembly begins. Since we have already applied diginorm to our reads, the result will be a negligible decrease in read counts prior to the assembly. We provide options to disable this digital normalization via the config file, but applying diginorm twice is not really a problem. For data sets with large numbers of reads, applying diginorm as a separate step as we have via khmer may decrease the memory requirements needed by the Trinity pipeline.

The ID for each transcript is output (version 2.2.0 to current) as follows, where the `TRINITY` is constant, the `DN2202` is an example of a variable contig/transcript ID, `c` stands for component, `g` gene and `i` isoform:

```
TRINITY_DN2202_c0_g1_i1
```

## Trinity Command

On the command line, the command eelpond runs is approximately:
```
Trinity --left left.fq \
  --right right.fq --seqType fq --max_memory 10G \
  --CPU 4
```

But we **highly recommend** you modify `max_memory` and `CPU` to fit your data and compute resources.

## Customizing Trinity parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a trinity configfile you can modify, run:
```
./run_eelpond trinity.yaml trinity --build_config
```
The output should be a small `yaml` configfile that contains:
```
####################  trinity  ####################
trinity:
  input_kmer_trimmed: true
  input_trimmomatic_trimmed: false
  add_single_to_paired: false # would you like to add the orphaned reads to the trinity assembly?
  max_memory: 30G
  seqtype: fq
  extra: ''
```
We recommend using kmer-trimmed reads as input. If both `input_kmer_trimmed` and `input_trimmomatic_trimmed` are `False`, we will just use raw reads from the `samples.tsv` file. 

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra trinity parameters, e.g.:
```
  extra: '--someflag someparam --someotherflag thatotherparam'
```
Or in Trinity params:
```
  extra: '--no_normalize_reads' # to turn off Trinity's digital normalization steps 
```
Override default params by modifying any of these lines, and placing them in the config file you're using to run `eelpond`. Here, we just generated params for `trinity`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file, e.g. `./run_eelpond my-workflow.yaml full --build_config` and edit parameters there.


## Trinity eelpond rule

We wrote a [Trinity snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/trinity.html) to run Trinity.

For snakemake afficionados, here's the basic structure of the trinity eelpond rule. Directories and parameters are specified via the configfile, for more, see the rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/trinity/trinity.rule).

```
rule trinity:
    input:
        unpack(get_assembly_input)
    output:
        fasta = join(assembly_dir,"trinity_out_dir/Trinity.fasta"),
        gene_trans_map = join(assembly_dir,"trinity_out_dir/Trinity.fasta.gene_trans_map"),
    message:
        """--- Assembling read data with Trinity --- """
    params:
        # optional parameters
        seqtype=assembly_params.get('seqtype', 'fq'),
        extra=assembly_params.get('extra', '')
    threads: 4
    log: join(LOGS_DIR, 'trinity/trinity.log')
    conda: "trinity-env.yaml"
	script: "trinity-wrapper.py"

rule rename_trinity_fasta:
    input: rules.trinity.output.fasta
    output: join(assembly_dir, BASE + assembly_extension + '.fasta')
    log: join(LOGS_DIR, 'trinity/cp_assembly.log')
    shell: ("cp {input} {output}") 

rule rename_trinity_gene_trans_map:
    input: rules.trinity.output.gene_trans_map
    output: join(assembly_dir, BASE + assembly_extension + '.fasta.gene_trans_map')
    log: join(LOGS_DIR, 'trinity/cp_gt_map.log')
    shell: ("cp {input} {output}") 
```

