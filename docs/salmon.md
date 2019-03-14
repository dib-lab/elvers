# Quantification with Salmon

We can use [Salmon](http://salmon.readthedocs.org/en/latest/) to quantify expression. Salmon is a (relatively) new breed of software for quantifying RNAseq reads that is both really fast and takes transcript length into consideration ([Patro et al. 2015](https://doi.org/10.1038/nmeth.4197)).

## Quickstart

We recommend that you run salmon quantification via the "default" [Eel Pond workflow](eel_pond_workflow.md) or the [quantify subworkflow](assemble.md). See "Advanced Usage" below for running salmon as a standalone rule.

## Salmon Commands 

There are two commands for salmon, `salmon index` and `salmon quant`. The first command, `salmon index` will index the transcriptome:

```
salmon index --index nema --transcripts nema_trinity.fasta --type quasi
```

And the second command, `salmon quant` will quantify the trimmed reads (not diginormed) using the transcriptome. For each pair of reads for a sample, we run:

```
salmon quant -i nema -l A -1 <(gunzip -c $R1) -2 <(gunzip -c $R2) -o ${sample_name}_quant
```

Both indexing the transcriptome and running quantification are integrated as rules in the elvers workflow, so the whole process happens in an automated fashion.

## Modifying Params for Salmon

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `salmon` rule, run
```
elvers config salmon --print_params
```
This will print the following:
```
  ####################  salmon  ####################
salmon:
  input_trimmomatic_trimmed: True
  index_params:
    extra: ''
  quant_params:
    libtype: A
    extra: '' 
  #####################################################
```
If you set `input_trimmomatic_trimmed: False` in the salmon parameters, then salmon will use your raw input data instead of trimming first. Using trimmed data as input is recommended, this is just if you're pre-trimmed with another program!

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra parameters.In salmon, both `index` and `quantification` steps can accept an `extra` param. See the [Salmon documentation](http://salmon.readthedocs.org/en/latest/) to learn more about the parameters you can pass into `salmon`.

Be sure the modified lines go into the config file you're using to run `elvers` (see [Understanding and Configuring Workflows](configure.md)).

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Salmon will output files in the `quant` subdirectory of this output directory. Each sample will have its own directory, and the two most interesting files will be the `salmon_quant.log` and `quant.sf` files. The former contains the log information from running salmon, and the latter contains the transcript count data.

A `quant.sf` file will look something like this.
```
Name                  Length    EffectiveLength    TPM    NumReads
TRINITY_DN2202_c0_g1_i1    210    39.818    2.683835    2.000000
TRINITY_DN2270_c0_g1_i1    213    41.064    0.000000    0.000000
TRINITY_DN2201_c0_g1_i1    266    69.681    0.766816    1.000000
TRINITY_DN2222_c0_g1_i1    243    55.794    2.873014    3.000000
TRINITY_DN2291_c0_g1_i1    245    56.916    0.000000    0.000000
TRINITY_DN2269_c0_g1_i1    294    89.251    0.000000    0.000000
TRINITY_DN2269_c1_g1_i1    246    57.479    0.000000    0.000000
TRINITY_DN2279_c0_g1_i1    426    207.443    0.000000    0.000000
TRINITY_DN2262_c0_g1_i1    500    280.803    0.190459    1.000912
TRINITY_DN2253_c0_g1_i1    1523    1303.116    0.164015    4.000000
TRINITY_DN2287_c0_g1_i1    467    247.962    0.000000    0.000000
TRINITY_DN2287_c1_g1_i1    325    113.826    0.469425    1.000000
TRINITY_DN2237_c0_g1_i1    306    98.441    0.542788    1.000000
TRINITY_DN2237_c0_g2_i1    307    99.229    0.000000    0.000000
TRINITY_DN2250_c0_g1_i1    368    151.832    0.000000    0.000000
TRINITY_DN2250_c1_g1_i1    271    72.988    0.000000    0.000000
TRINITY_DN2208_c0_g1_i1    379    162.080    1.978014    6.000000
TRINITY_DN2277_c0_g1_i1    269    71.657    0.745677    1.000000
TRINITY_DN2231_c0_g1_i1    209    39.409    0.000000    0.000000
TRINITY_DN2231_c1_g1_i1    334    121.411    0.000000    0.000000
TRINITY_DN2204_c0_g1_i1    287    84.121    0.000000    0.000000
```

## More on Salmon
For further reading, on salmon see

  * Intro blog post: http://robpatro.com/blog/?p=248
  * A 2016 blog post evaluating and comparing methods [here](https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/)
  * Salmon github repo [here](https://github.com/COMBINE-lab/salmon)
  * https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/salmon.rst
  * http://angus.readthedocs.io/en/2016/rob_quant/tut.html
  * https://2016-aug-nonmodel-rnaseq.readthedocs.io/en/latest/quantification.html

## Advanced Usage: Running Salmon as a standalone rule

You can run salmon as a standalone rule, instead of withing a larger `elvers` workflow. However, to do this, you need to make sure the input files are available.

For salmon, you need both 1) an assembly, and 2) trimmed input files. The assembly can be generated via another workflow, or passed to `elvers` via the configfile.

Specifying an assembly:
    1) If you've alread run read trimming and want to use a Trinity assembly generated via `elvers`, you can run: 
    ```
    elvers my_config assemble salmon
    ```
    If you've already run the assembly, `elvers` will just use this info to locate that assembly.

    2) Alternatively, you can input an assembly via the [assemblyinput](assemblyinput.md) utility rule:
    ```
    elvers assemblyinput salmon
     ```
    with an assembly in your `yaml` configfile, e.g.:
    ```
    assemblyinput:
      assembly: examples/nema.assembly.fasta
      gene_trans_map:  examples/nema.assembly.fasta.gene_trans_map #optional
      assembly_extension: '_input'
    ```
    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have a gene to transcript map, please specify it as well. If not, delete this line from  your `config`. The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Feel free to use `_input`, as specified above, or pick something equally simple yet more informative. **Note: Please don't use additional underscores (`_`) in this extension!**. For more details, see the [assemblyinput documentation](assemblyinput.md).

Specifying input reads:

    If you haven't yet run read trimming, you'll also need to run those steps:
    ```
    elvers myconfig get_data trimmomatic salmon
    ```
    Or if you have set `input_trimmomatic_trimmed: False`:
    ```
    elvers myconfig get_data salmon
    ```

## Snakemake Rule

We wrote snakemake wrappers to run [salmon index](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/salmon/index.html) and [salmon quant](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/salmon/quant.html).

For snakemake afficionados, see the Salmon rule on [github](https://github.com/dib-lab/elvers/blob/master/rules/salmon/salmon.rule).
