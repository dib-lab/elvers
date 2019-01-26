# Understanding and Configuring Eelpond Workflows

`eelpond` is designed to facilitate running standard workflows and analyses for sequence data. It integrates snakemake rules for commonly used tools, and provides several end-to-end protocols for analyzing RNAseq data. Workflows are highly customizable, and all command-line options for each tool are available for modification via modification of the configuration file.

## Specify Input Data

To run any workflow, we need to get reads (and any assemblies already generated) into the right format for `eelpond`. We can either (1) link the data from another location on your machine, or (2) download the data via http or ftp. This is done through a utility rule called `get_data`.

To specify input data, we generate a `my-samples.tsv` file, with either file paths or download links

If we had our test data within the `data` folder in the main `eelpond` directory, the file would look like this:

    ```
    sample  unit    fq1 fq2 condition
    0Hour   001 data/0Hour_ATCACG_L002_R1_001.extract.fastq.gz data/0Hour_ATCACG_L002_R2_001.extract.fastq.gz time0
    0Hour   002 data/0Hour_ATCACG_L002_R1_002.extract.fastq.gz data/0Hour_ATCACG_L002_R2_002.extract.fastq.gz time0
    6Hour   001 data/6Hour_CGATGT_L002_R1_001.extract.fastq.gz data/6Hour_CGATGT_L002_R2_001.extract.fastq.gz time6
    6Hour   002 data/6Hour_CGATGT_L002_R1_002.extract.fastq.gz data/6Hour_CGATGT_L002_R2_002.extract.fastq.gz time6
    ```

If we want to download the test data instead:

    ```
    sample  unit    fq1 fq2 condition
    0Hour   001 https://osf.io/vw4dt/download   https://osf.io/b47s2/download   time0
    0Hour   002 https://osf.io/92jr6/download   https://osf.io/qzea8/download   time0
    6Hour   001 https://osf.io/ft3am/download   https://osf.io/jqmsx/download   time6
    6Hour   002 https://osf.io/rnkq3/download   https://osf.io/bt9vh/download   time6
    ```

**Note that proper formatting means all columns must be separated by tabs, not spaces!**

If you'd like to start from a working version, copy the sample data:

```
cp examples/nema.samples.tsv my_samples.tsv
```

Note: at the moment, `eelpond` assumes **all input data is gzipped**.

## Choosing and running a workflow

The default workflow is the [Eel Pond RNAseq workflow](eel_pond_workflow.md), but we offer several other end-to-end workflows, as well as a number of "subworkflows" that facilitate running (or re-running) certain steps of the workflows. 

To see the available workflows and subworkflows, run:
```
./run_eelpond examples/nema.yaml -w 
```

## Configuring Programs in the workflow

For any workflow, we need to provide a configuration file that specifies the path to your 

To get a `preprocess` config you can modify, run:

```
./run_eelpond my_workflow.yaml --build_config
```

The output should be a `yaml` configfile. At the top, you should see:

```
  ####################  Eelpond Pipeline Configfile  ####################
basename: eelpond
experiment: _experiment1
samples: samples.tsv
```

## Customizing the Output Location and Names

First, change the `samples.tsv` name to your `my_samples.tsv` file. Then, modify the basename and any experiment info you'd like to add. The default output directory will be: `basename_experiment_out` within the main `eelpond` directory. If you'd like, you can add one more parameter to the top section: `out_path: OUTPUT_PATH`, if you'd like the output to go somewhere other than the `eelpond` directory. Note the basename and experiment are still used to determine the output directory name.


## Customizing program parameters:

Below this section, you should see some parameters for each program run in this workflow. For example, here's the first few programs: a utility to download or link your data, quality trim with trimmomtic, and assess quality with fastqc. Program parameters do not always show up in order in this file - order in this file does not affect program run order.

```
get_data:
  download_data: false
  use_ftp: false
trimmomatic:
  adapter_file:
    pe_path: ep_utils/TruSeq3-PE-2.fa
    se_path: ep_utils/TruSeq3-SE.fa
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25
  extra: ''
fastqc:
  extra: ''
```

Override default params for any program by modifying the values for any of the paramters under that program name. For example, if you'd like to download data instead of link it from another location on your machine, modify `download_data` to `True` under the `get_data` program. We provide an `extra` parameter wherever possible to give you access to additional command-line parameters for each program that we have not specifically enabled changes for.

For example, under the `trimmomatic` section, you can modify the "extra" param to pass any extra trimmomatic parameters, e.g.:

```
trimmomatic:
  extra: 'HEADCROP:5' # to remove the first 5 bases at the front of the read.
```

For more on what parameters are available, see the docs for each specific program or utility rule.

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md)
  - [khmer](khmer.md)
  - [trinity](trinity.md)
  - [dammit](dammit.md)
  - [salmon](salmon.md)
  - [sourmash](sourmash.md)
