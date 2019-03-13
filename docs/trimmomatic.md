# Trimmomatic

We use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.36) to trim off residual Illumina adapters that were left behind after demultiplexing.

Here we use a set TruSeq Illumina adapters. However, if running this on your own data and you know you have different adapters, you'll want to input them in the configfile (see `params` section, below). If you're using the right adapters,  you should see that some of the reads are trimmed; if they’re not, you won’t see anything get trimmed.

See this excellent paper on trimming parameters by [MacManes 2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full).

## Trimmomatic Command

Based on recommendations from [MacManes 2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full), we use this command by default:

```
TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:15 \
    MINLEN:25
```

However, the trimming command can be extensively modified via the configfile. Here's how the parameters for the command above look in our config:

```
trimmomatic:
  trim_cmd:  "ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25"
  extra: ''
```
## Quickstart

Run Trimmomatic via the "default" [Eel Pond workflow](eel_pond_workflow.md) or via the [preprocess subworkflow](preprocess.md). To run Trimmomatic as a standalone program, see "Advanced Usage" section below.

## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Trimmomatic will output files in the `preprocess` subdirectory of this output directory. All outputs will contain `*.trim.fq.gz`.

## Modifying Params for Trimmomatic:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `trimmomatic` rule, run
```
./run_elvers config trimmomatic --print_params
```

In here, you'll see a section for "trimmomatic" parameters that looks like this:

```
  ####################  trimmomatic  ####################
trimmomatic:
  adapter_file:
    pe_path: ep_utils/TruSeq3-PE-2.fa
    se_path: ep_utils/TruSeq3-SE.fa
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
  extra: ''
  #####################################################
```

Override default params by modifying these lines. In addition to changing parameters we've specifically enabled, you can modify the "extra" param to pass any extra trimmomatic parameters, e.g.:

```
  extra: 'HEADCROP:5' # to remove the first 5 bases at the front of the read.
```
See the [Trimmomatic documentation](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for parameters and options you can pass to Trimmomatic. Be sure the modified lines go into the config file you're using to run `elvers` (see [Understanding and Configuring Workflows](configure.md)).


## Advanced Usage: Running Trimmomatic as a standalone rule

You can run trimmomatic as a standalone rule, instead of withing a larger `elvers` workflow. However, to do this, you need to make sure the input files are available.

For Trimmomatic, the input files are your input data - either downloaded or linked into the `input_data` directory via `get_data`.

If you've already done this, you can run:
```
./run_elvers my_config trimmomatic
```
If not, you can run both at once to make sure trimmomatic can find its input files:
```
./run_elvers my_config get_data trimmomatic
```


## Snakemake Rule 

We use a local copies of the [trimmomatic snakemake wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/trimmomatic.html) to run Trimmomatic.

For snakemake afficionados, see the trimmomatic rules on [github](https://github.com/dib-lab/elvers/blob/master/rules/trimmomatic/trimmomatic.rule).
