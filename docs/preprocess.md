# Preprocess

The "preprocess" workflow is a subworkflow that conducts read quality trimming.

# Quickstart

To run the preprocess workflow, run: 

```
./run_eelpond nema-test preprocess
```

At the moment, this workflow consists of:
 
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md), run on both pre-trim and post-trim data


## Specify Input Data

To run any workflow, we need to get reads (and any assemblies already generated) into the right format for `eelpond`. We can either (1) link the data from another location on your machine, or (2) download the data via http or ftp. This is done through a utility rule called `get_data`. 

To specify input data, we generate a `my-samples.tsv` file, with either file paths or download links

If we had our test data within the `nema_testdata` folder in the main `eelpond` directory, the file would look like this:

    ```
    sample  unit    fq1 fq2 condition
    0Hour   001 nema_testdata/0Hour_ATCACG_L002_R1_001.extract.fastq.gz nema_testdata/0Hour_ATCACG_L002_R2_001.extract.fastq.gz time0
    0Hour   002 nema_testdata/0Hour_ATCACG_L002_R1_002.extract.fastq.gz nema_testdata/0Hour_ATCACG_L002_R2_002.extract.fastq.gz time0
    6Hour   001 nema_testdata/6Hour_CGATGT_L002_R1_001.extract.fastq.gz nema_testdata/6Hour_CGATGT_L002_R2_001.extract.fastq.gz time6
    6Hour   002 nema_testdata/6Hour_CGATGT_L002_R1_002.extract.fastq.gz nema_testdata/6Hour_CGATGT_L002_R2_002.extract.fastq.gz time6
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
cp nema_samples.tsv my_samples.tsv
```

Note: at the moment, `eelpond` assumes **all input data is gzipped**.

## Building a Configuration File 

To get a `preprocess` config you can modify, run:

```
./run_eelpond preprocess.yaml preprocess --build_config
```

The output should be a small `yaml` configfile that contains:

```

  ####################  Eelpond Pipeline Configfile  ####################
basename: eelpond
experiment: _experiment1
samples: samples.tsv

  ####################  preprocess  ####################
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

## Customizing the Output Location and Names

First, change the `samples.tsv` name to your `my_samples.tsv` file. Then, modify the basename and any experiment info you'd like to add. The default output directory will be: `basename_experiment_out` within the main `eelpond` directory. If you'd like, you can add one more parameter to the top section: `out_path: OUTPUT_PATH`, if you'd like the output to go somewhere other than the `eelpond` directory. Note the basename and experiment are still used to determine the output directory name.


## Customizing program parameters:

Override default params for any program by modifying the values for any of the paramters under that program name. For example, if you'd like to download data instead of link it from another location on your machine, modify `download_data` to `True` under the `get_data` program. We provide an `extra` parameter wherever possible to give you access to additional command-line parameters for each program that we have not specifically enabled changes for.

For example, under the `trimmomatic` section, you can modify the "extra" param to pass any extra trimmomatic parameters, e.g.:

```
trimmomatic:
  extra: 'HEADCROP:5' # to remove the first 5 bases at the front of the read.
```

For more on what parameters are available, see the docs for each specific program or utility rule:

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md)

Here, we just generated params for `preprocess`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file, e.g. `./run_eelpond my-workflow.yaml eel_pond --build_config` and edit parameters there. If you forget, you can always generate an additional configfile and copy in the program-specific parameters that you need.

