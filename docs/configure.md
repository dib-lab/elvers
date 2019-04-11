# Understanding and Configuring Eelpond Workflows

`elvers` is designed to facilitate running standard workflows and analyses for sequence data. It integrates snakemake rules for commonly used tools, and provides several end-to-end protocols for analyzing RNAseq data. Workflows are highly customizable, and all command-line options for each tool are available for modification via the configuration file.


**`elvers`** uses the `yaml` _Yet Another Markup Language_ format to specify data inputs and modify run parameters. The only required information for this file is the location of the data inputs (reads, assembly, or both). 

## Specifying Input Data

Let's start with the input data. There are two types of data that can go into elvers: read data (gzipped fastq files), and assembly files (fasta) files. 

### Read Input

To specify input data, we need to build a tab-separated samples file, e.g. `my-samples.tsv`. This file tells `elvers` a name for each samples, and provides a location for the fastq files (local file path, or downloadable link).

If we had our test data within the `data` folder in the main `elvers` directory, the file would look like this:

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

Now, you need to provide the name and location of this file to `elvers`. We do this in the `config.yaml` file, like so:
```
get_data:
  samples: /path/to/my-samples.tsv
```

The default functionality is to link data from another location on your computer into `elvers`'s output directory. However, if you want to download data instead, you'll need to provide links in the `my-samples.tsv` file, and add a few lines to your `config.yaml`:

```
get_data:
  samples: /path/to/my-samples.tsv
  download_data: True 
  use_ftp: False # set to true if you want to use FTP instead of HTTP
```

**About the TSV and Input Reads:**
  
  - The sample names provided here are used to name file outputs throughout the workflow. 
  - The "unit" column is designed to facilitate combining samples that were sequenced over multiple lanes or in different batches. If you do not have "unit" information, please add a short placeholder, such as "a".
  - for now, the column headers must be: `sample`, `unit`, `fq1`, `fq2`, `condition`. Additional headers are not a problem, but will not be used. 
  - If you have single-end samples, please be sure to include the `fq2` column, but leave the column blank
  - the `condition` column is used for differential expression comparisions in DESeq2. The "condition" values can be anything, but they need to correspond with the "contrast" information you pass into DESeq2. If using the `diffexp` workflow, see our docs [here](diffexp.md) and be sure to add the differential expression information into your `yaml` configuration file.
  - Formatting the `tsv` can be a bit annoying. It's slightly easier if you start from a working version by copying the sample data to a new file (see below).
  - At the moment, `elvers` assumes **all input data is gzipped**, so please input gzipped data!


If you'd like to start from a working version, copy (and then modify) the sample data:

```
cp examples/nema.samples.tsv my_samples.tsv
```

To run any read-based workflow, we need to get reads (and any assemblies already generated) into the right format for `elvers`. We can either (1) link the data from another location on your machine, or (2) download the data via http or ftp. This is done through a utility rule called `get_data`.

### Reference Input

If you're starting a new `elvers` run using a reference file (even if you have previously built a _de novo_ assembly via `elvers`), you need to help `elvers` find that assembly.


**Scenario 1: You're starting from your own reference file:**

You need to provide the reference in your `config.yaml` file:
```
get_reference:
  reference: input reference fasta REQUIRED
  gene_trans_map: OPTIONAL: provide a gene to transcript map for input transcriptome
  reference_extension: '_input' OPTIONAL, changes naming
  download_ref: the reference entry above is a link that needs to be downloaded
  use_ftp: download via ftp instead of http
```
Once you add this to your configfile (e.g. `my_config.yaml`), you can run a reference/assembly-based workflow such as annotate.
```
elvers my_config.yaml annotate
```
For more details on reference specification, see the [get_reference documentation](get_reference.md). For annotation configuration, see below.

**Scenario 2: You've previously run an assembly program via `elvers`:**

You _can_ provide the built reference in the same manner as above. However, if you're running more workflows in the same directory, you can also just specify the name of the assembly program that you used to generate the assembly. This will *not* rerun the assembly (unless you provide new input files). Instead, this will allow `elvers` to know where to look for your previously-generated reference     file. Because we enable multiple referene generation programs, we don't want to assume
which reference you'd like to use for downstream steps (in fact, if you provide multiple references, `elvers` will run the downstream steps on all references, assuming you provide unique               `reference_extension` parameters so that the references are uniquely named.

Example: You've previously run the trinity assembly, and want to annotate it.

```
elvers examples/nema.yaml assemble annotate
```
Here, the `assemble` workflow just enables `elvers` to locate your assembly file for `annotate`.


## Choosing and running a workflow

We offer a number of workflows, including end-to-end workflows that conduct assembly through differential expression analysis. The default workflow is the [Eel Pond RNAseq workflow](eel_pond_workflow.md), which conducts *de novo* transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command. 

### Available Workflows

Currently, all workflows require a properly-formatted read inputs `tsv` file as input. Some workflows, e.g. `annotation`, can work on either on a  _de novo_ transcriptome generated by `elvers`, or on previously-generated assemblies. To add an assembly as input, specify it via `get_reference` in the `yaml` config file, as described above. 


**workflows**

  - preprocess: Read Quality Trimming and Filtering (fastqc, trimmomatic)
  - kmer_trim: Kmer Trimming and/or Digital Normalization (khmer)
  - assemble: Transcriptome Assembly (trinity)
  - get_reference: Specify assembly for downstream steps
  - annotate : Annotate the transcriptome (dammit)
  - sourmash_compute: Build sourmash signatures for the reads and assembly (sourmash)
  - quantify: Quantify transcripts (salmon)
  - diffexp: Conduct differential expression (DESeq2)
  - plass_assemble: assemble at the protein level with PLASS
  - paladin_map: map to a protein assembly using paladin

**end-to-end workflows:**

  - **default**: preprocess, kmer_trim, assemble, annotate, quantify
  - **protein assembly**: preprocess, kmer_trim, plass_assemble, paladin_map


You can see the available workflows (and which programs they run) by using the `--print_workflows` flag:
```
elvers examples/nema.yaml --print_workflows
```

Each included tool can also be run independently, if appropriate input files are provided. This is not always intuitive, so please see our documentation for running each tools for details (described as "Advanced Usage"). To see all available tools, run:

```
elvers examples/nema.yaml --print_rules
```

### Configuring Parameters for a workflow

For any workflow, we need to provide a configuration file that specifies the path to your samples file or get_reference (discussed above).

We can generate this file either by:

  1. Adding just the desired parameters
  2. Allowing `elvers` to build a (long) full configfile for us, and modifying as desired.


**Option 1: Adding just the required parameters**

The configuration file primarily provides the location of the input data and/or input assembly. The simplest config file contains just this information.

```
get_data:
  samples: samples.tsv
```
or 
```
get_reference:
  reference: assembly.fasta
```


There are a few other options we can add to customize the name of the output directory and files.

  - `basename: NAME`: helps determine file names and output directory (by default: `BASENAME_out`)
  - `experiment: EXPERIMENT`: some additional "experiment" info to add to the output directory name ( outdir: `BASENAME_EXPERIMENT_out`)
  - `out_path: /full/path`: if you want to redirect the output to some location *not* under the `elvers` directory.
  - Finally, to specify a specific set of workflows or tools to use, add `workflows` to your `yaml` file. In this case, we just want to run `fastqc` and `trimmomatic`:

```
workflows: 
  - fastqc
  - trimmomatic
``` 


Now, if you'd like to run any particular program with non-default parameters, or you're running differential expression analysis, you'll need to add some info to the config. For any (each) program, follow this format to see program params:

```
elvers config.yaml <PROGRAM> --print_params
```

For example, for deseq2:
```
elvers config.yaml deseq2 --print_params
```

Then copy and paste the parameters that show up in your terminal into your config file. Please see each program's documentation for additional info on what (and how) to modify each program. For example, if running an assembly, we definitely recommend modifying the `max_memory` parameter of Trinity.


**Option 2: Use `elvers` to add all parameters for our workflow of choice**

To get a configfile for the default "eel pond" workflow, that you can modify, run:

```
elvers my_workflow.yaml --build_config
```

The output should be a `yaml` configfile. At the top, you should see:

```
#
# elvers workflow configuration
#
basename: elvers
experiment: _experiment1
```

First, change the `samples.tsv` name to your `my_samples.tsv` file. Then, modify the basename and any experiment info you'd like to add. The default output directory will be: `basename_experiment_out` within the main `elvers` directory. If you'd like, you can add one more parameter to the top section: `out_path: OUTPUT_PATH`, if you'd like the output to go somewhere other than the `elvers` directory. Note the basename and experiment are still used to determine the output directory name.


**Customizing program parameters:**

Below this section, you should see some parameters for each program run in this workflow. For example, here's the first few programs: a utility to download or link your data, quality trim with trimmomtic, and assess quality with fastqc. Program parameters do not always show up in order in this file - order in this file does not affect program run order.

```
get_data:
  samples: samples.tsv
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

For more on what parameters are available, see the docs for each specific program or utility rule under the "Available Workflows" --> "Programs Used" navigation tab.
