# Protein Assembly with PLASS

The basic idea with any transcriptome or metatranscriptome assembly is you feed in your reads and you get out a bunch of *contigs* that represent transcripts, or stretches of RNA present in the reads. You run a transcriptome assembly program using the adapter & k-mer trimmed reads as input and get out a pile of assembled RNA. 

Most assemblers work in *nucleotide* space. This means that they find direct overlaps between the As, Ts, Cs, and Gs in the reads and use information about those overlaps to make contigs. Although this is a powerful method, it often fails when:

  -  There is sequencing errors
  -  There are repetitive regions
  -  There is strain variation

When assembly fails, contigs either break or don't assemble at all. 

Given that nucleotide sequences are far more variable than protein sequences (third base pair wobble in strain variation in particular), assembling in amino acid space can overcome a lot of the difficulties encountered when assembling in nucleotide space. 

[PLASS](https://plass.mmseqs.org) is a new assembler that assembles in amino acid space. Unlike many other assemblers, including [Trinity](trinity.md), it does not have built in error correction and so it is best to adapter and k-mer trim the reads before using it. It sometimes performs better than nucleotide assemblers and so is good to test on samples to see what type of improvement it can give. 

## Quickstart

Run PLASS via the [Protein Assembly workflow](protein_assembly_workflow.md) or via the [plass_assemble subworkflow](plass_assemble.md). These workflows will run preprocessing and kmer-trimming for you prior to assembly, and may run additional downstream steps. To run PLASS as a standalone program, see "Advanced Usage" section below.
```
./run_eelpond examples/nema.yaml plass_assemble
```

## PLASS Command

On the command line, the command eelpond runs is approximately:
```
plass assemble input_1.fq input_2.fq \
  outputi_plass.fasta --threads snakemake.threads extra_params
```

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

PLASS will output files in the `assembly` subdirectory of this output directory. The fasta file will be `BASENAME_plass.fasta`.

## Modifying Params for PLASS:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `plass` rule, run
```
./run_eelpond config plass --print_params
```
This will print the following:
```
  ####################  plass  ####################
plass:
  input_kmer_trimmed: true
  input_trimmomatic_trimmed: false
  add_single_to_paired: false # would you like to add the orphaned reads to the plass assembly?
  extra: '' 
  #####################################################
```
Within the [Protein Assembly workflow](protein_assembly_workflow.md) or the [plass_assemble subworkflow](plass_assemble.md), these options enable you to choose kmer-trimmed, quality-trimmed, or raw sequencing data as input. We recommend using kmer-trimmed reads as input. If both `input_kmer_trimmed` and `input_trimmomatic_trimmed` are `False`, we will just use raw reads from the `samples.tsv` file.

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra plass parameters,  e.g.:
```
  extra: '--someflag someparam --someotherflag thatotherparam'
```
See the [PLASS documentation](https://plass.mmseqs.org) to learn more about the parameters you can pass to `PLASS`.

Be sure the modified lines go into the config file you're using to run `eelpond` (see [Understanding and Configuring Workflows](configure.md)).

## Advanced Usage: Running PLASS as a standalone rule

You can run plass as a standalone rule, instead of withing a larger `eelpond` workflow. However, to do this, you need to make sure the input files are available.

For plass, the default input files are kmer-trimmed input data (e.g. output of khmer).

If you've already done this, you can run:
```
./run_eelpond my_config plass
```
If not, you can run the prior steps at the same time to make sure khmer can find these input files:
```
./run_eelpond my_config get_data trimmomatic khmer plass
```

## Snakemake Rules

We wrote a [PLASS snakemake wrapper](https://github.com/dib-lab/eelpond/blob/master/rules/plass/plass-wrapper.py) to run PLASS via snakemake. This wrapper has not yet been submitted to the snakemake-wrappers repository, but feel free to use it as needed.

For snakemake afficionados, see our PLASS rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/plass/plass.rule).
