# Assembling with PLASS


## Quickstart: running PLASS via eelpond:

```
./run_eelpond nema-test plass_assemble
```
This will run preprocessing and kmer-trimming for you prior to assembly. 

## Protein Assembly with PLASS

The basic idea with any transcriptome or metatranscriptome assembly is you feed in your reads and you get out a bunch of *contigs* that represent transcripts, or stretches of RNA present in the reads. You run a transcriptome assembly program using the adapter & k-mer trimmed reads as input and get out a pile of assembled RNA. 

Most assemblers work in *nucleotide* space. This means that they find direct overlaps between the As, Ts, Cs, and Gs in the reads and use information about those overlaps to make contigs. Although this is a powerful method, it often fails when:

  -  There is sequencing errors
  -  There are repetitive regions
  -  There is strain variation

When assembly fails, contigs either break or don't assemble at all. 

Given that nucleotide sequences are far more variable than protein sequences (third base pair wobble in strain variation in particular), assembling in amino acid space can overcome a lot of the difficulties encountered when assembling in nucleotide space. 

[PLASS](https://plass.mmseqs.org) is a new assembler that assembles in amino acid space. Unlike many other assemblers, including [Trinity](trinity.md), it does not have built in error correction and so it is best to adapter and k-mer trim the reads before using it. It sometimes performs better than nucleotide assemblers and so is good to test on samples to see what type of improvement it can give. 

## PLASS Command

On the command line, the command eelpond runs is approximately:
```
plass assemble input_1.fq input_2.fq \
  output.fasta --threads snakemake.threads extra_params
```

## Customizing PLASS parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a PLASS configfile you can modify, run:
```
./run_eelpond plass.yaml plass --build_config
```
The output should be a small `yaml` configfile that contains:
```
  ####################  plass  ####################
plass:
  input_kmer_trimmed: true
  input_trimmomatic_trimmed: false
  add_single_to_paired: false # would you like to add the orphaned reads to the plass assembly?
  extra: ''
```
We recommend using kmer-trimmed reads as input. If both `input_kmer_trimmed` and `input_trimmomatic_trimmed` are `False`, we will just use raw reads from the `samples.tsv` file.

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra plass parameters,  e.g.:
```
  extra: '--someflag someparam --someotherflag thatotherparam'
```
Override default params by modifying any of these lines, and placing them in the config file you're using to run `eelpond`. Here, we just generated params for `plass`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file, e.g. `./run_eelpond my-workflow.yaml full --build_config` and edit parameters there. 


## PLASS eelpond rule

We wrote a [PLASS snakemake wrapper](https://github.com/dib-lab/eelpond/blob/master/rules/plass/plass-wrapper.py) to run PLASS via snakemake. This wrapper has not yet been submitted to the snakemake-wrappers repository, but feel free to use it as needed.

For snakemake afficionados, see our PLASS rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/plass/plass.rule).





