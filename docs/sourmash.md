# Sourmash

sourmash is a command-line tool and Python library for computing MinHash sketches from DNA sequences, comparing them to each other, and plotting the results. This allows you to estimate sequence similarity between even very large data sets quickly and accurately. Please see the mash software and the mash paper [(Ondov et al., 2016)](http://biorxiv.org/content/early/2015/10/26/029827) for background information on how and why MinHash sketches work.

Sourmash is [dib-lab](http://ivory.idyll.org/lab/) software! Please see the [sourmash documentation](https://sourmash.readthedocs.io/en/latest/index.html) for more on sourmash. Sourmash 2.0 is coming soon. In the meantime, please cite [Brown and Irber, 2016](https://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9)

**At the moment we have only enabled _sourmash compute_ functionality.**

## Quickstart

Run Sourmash as part of the "default" [Eel Pond workflow](eel_pond_workflow.md) or via the [sourmash_compute subworkflow](sourmash_compute.md). At the moment, sourmash compute requires both an assembly and a set of reads as input. Please see the [sourmash_compute subworkflow](sourmash_compute.md) for how to run sourmash compute properly. 

## Sourmash Command

On the command line, the command eelpond runs for each file is approximately:
```
sourmash compute --scaled 1000 \
  -k 31 input_file -o output.sig
```

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Sourmash will output files in the `sourmash` subdirectory of this output directory. Sourmash signatures will have the same name as the file they're generated from, but end with `.sig` instead of `.fasta` or `.fq.gz`.


## Modifying Params for Sourmash:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `sourmash` rule, run
```
./run_eelpond config sourmash --print_params
```
This will print the following:
```
  ####################  sourmash  ####################
sourmash:
  k_size: 31
  scaled: 1000
  extra: '' 
  #####################################################
```
In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra sourmash parameters,  e.g.:
```
  extra: ' --track-abundance '
```
Be sure the modified lines go into the config file you're using to run `eelpond` (see [Understanding and Configuring Workflows](configure.md)). 

See the [sourmash documentation](https://sourmash.readthedocs.io/en/latest/index.html) to learn more about the parameters you can use with sourmash compute. 


## Sourmash eelpond rule

We use a slightly modified version of the [sourmash snakemake wrapper](https://github.com/dib-lab/eelpond/blob/master/rules/sourmash/sourmash-wrapper.py) to run Sourmash compute via snakemake. 

For snakemake afficionados, see our sourmash rules on [github](https://github.com/dib-lab/eelpond/blob/master/rules/sourmash/sourmash.rule).
