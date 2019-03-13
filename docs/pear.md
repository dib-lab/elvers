# Merging Read Pairs with PEAR

In order to map PE nucleotide reads to protein assemblies with `Paladin`, which cannot yet use paired end reads, we use `PEAR` to merge PE reads prior to mapping.

From the `PEAR` [documentation](https://cme.h-its.org/exelixis/web/software/pear/doc.html):

PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger. It is fully parallelized and can run with as low as just a few kilobytes of memory.

PEAR evaluates all possible paired-end read overlaps and without requiring the target fragment size as input. In addition, it implements a statistical test for minimizing false-positive results. Together with a highly optimized implementation, it can merge millions of paired end reads within a couple of minutes on a standard desktop computer.

## Quickstart: Running PEAR with elvers

We recommend you run `pear` as part of the [paladin_map](paladin_map.md) subworkflow

```
./run_elvers examples/nema.yaml paladin_map
```
This will run trimmomatic trimming prior to PEAR merging of paired end reads and then paladin mapping. If you'd like to just run `PEAR`, see "Advanced Usage" below.

## PEAR Command

On the command line, the command elvers runs for each set of files is approximately:
```
pear -f forward_reads -r reverse_reads \ 
  -p p_value -j snakemake.threads -y max_memory \  
  extra -o output_basename 
```

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

PEAR will output files in the `preprocess/pear` subdirectory of this output directory. PEAR output will have the same sample name as input, but end with `.pear.fq.gz`. 

## Modifying Params for PEAR:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `PEAR` rule, run
```
./run_elvers config pear --print_params
```
This will print the following:
```
  ####################  pear  ####################
pear:
  input_kmer_trimmed: false
  input_trimmomatic_trimmed: true
  max_memory: 4G
  pval: 0.01
  extra: ''
  #####################################################
```
**Please modify the `max_memory` parameter to fit the needs of your reads and system.**
Within the [Protein Assembly workflow](protein_assembly_workflow.md) or the [paladin_map subworkflow](paladin_map.md), these options enable you to choose kmer-trimmed, quality-trimmed, or raw sequencing data as input. We recommend using quality-trimmed reads as input. If both `input_kmer_trimmed` and `input_trimmomatic_trimmed` are `False`, we will just use raw reads from the `samples.tsv` file.

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass in additional parameters,  e.g.:

```
  extra: ' --some_param that_param '
```
Please see the [PEAR documentation](https://cme.h-its.org/exelixis/web/software/pear/doc.html) for info on the params you can pass into PEAR.

Be sure the modified lines go into the config file you're using to run `elvers` (see [Understanding and Configuring Workflows](configure.md)).


## Advanced Usage: Running PEAR as a standalone rule

You can run pear as a standalone rule, instead of withing a larger `elvers` workflow. However, to do this, you need to make sure the input files are available.

For pear, the default input files are quality-trimmed input data (e.g. output of trimmomatic).

If you've already done this, you can run:
```
./run_elvers my_config pear
```
If not, you can run the prior steps at the same time to make sure pear can find these input files:
```
./run_elvers my_config get_data trimmomatic pear
```
Note, `PEAR` only works on paired end reads, as there's nothing to do for single end reads!

## Snakemake Rule

We wrote a new [PEAR snakemake wrapper](https://github.com/dib-lab/elvers/blob/master/rules/pear/pear-wrapper.py) to run PEAR via snakemake. This wrapper has not been added to the official snakemake-wrappers repo yet, but feel free to use as needed.

For snakemake afficionados, see our pear rule on [github](https://github.com/dib-lab/elvers/blob/master/rules/pear/pear.rule).
