# Khmer k-mer trimming and (optional) diginorm

Before running transcriptome assembly, we recommend doing some kmer spectral error trimming on your dataset, and if you have lots of reads, also performing digital normalization. This has the effect of reducing the computational cost of assembly without negatively affecting the quality of the assembly. We use [khmer](https://khmer.readthedocs.io/) for both of these tasks.

# Khmer Command

Here's the command as it would look on the command line:
```
(interleave-reads.py sample_1.trim.fq sample_2.trim.fq )| \\
(trim-low-abund.py -V -k 20 -Z 18 -C 2 - -o - -M 4e9 --diginorm --diginorm-coverage=20) | \\
(extract-paired-reads.py --gzip -p sample.paired.gz -s sample.single.gz) > /dev/null
```

Trimmed reads are used as input. `trim-low-abund.py ` trims low-abundance reads to a coverage of 18. Here, we also perform digital normalization to a *k*-mer (*k* = 20) coverage of 20.The output files are the remaining reads, grouped as pairs and singles (orphans). For more on `trim-low-abund`, see this [recipe](https://khmer-recipes.readthedocs.io/en/latest/007-variable-coverage-trimming/),

Finally, since Trinity expects separate `left` and `right` files, we use `split-paired-reads.py` to split the interleaved pairs into two files.
```
split-paired-reads.py sample.paired.gz
```

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Khmer will output quality control files in the `preprocess` subdirectory of this output directory. All outputs will contain `*.khmer.fq.gz`.


## Modifying Params for Khmer:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](about_and_configure.md)).

To see the available parameters for the `trimmomatic` rule, run
```
./run_eelpond config trimmomatic --print_params
```

In here, you'll see a section for "trimmomatic" parameters that looks like this:

```
    ####################  khmer  ####################
khmer:
  C: 3
  Z: 18
  coverage: 20
  diginorm: true
  ksize: 20
  memory: 4e9
  #####################################################
```

See the [Khmer documentation]([khmer](https://khmer.readthedocs.io/)) to learn more about these parameters. Be sure the modified lines go into the config file you're using to run `eelpond` (see [Understanding and Configuring Workflows](about_and_configure.md)).

## Advanced Usage: Running Khmer as a standalone rule

You can run khmer as a standalone rule, instead of withing a larger `eelpond` workflow. However, to do this, you need to make sure the input files are available.

For khmer, the input files are trimmed input data (e.g. output of trimmomatic). 

If you've already done this, you can run:
```
./run_eelpond my_config khmer
```
If not, you can run the prior steps at the same time to make sure khmer can find these input files: 
```
./run_eelpond my_config get_data trimmomatic khmer
```


## Snakemake Rule

For snakemake afficionados, see the khmer rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/khmer/khmer.rule).
