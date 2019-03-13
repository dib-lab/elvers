# Preprocess Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting withing `elvers`. The "preprocess" subworkflow conducts read quality trimming. At the moment, this workflow consists of:
 
  - [get_data](get_data.md) - an `elvers` utility
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md), run on both pre-trim and post-trim data


# Quickstart

To run the preprocess subworkflow, run: 

```
./run_elvers examples/nema.yaml preprocess
```

## Configuring the preprocess workflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](configure.md).

If you want to add the `preprocess` program parameters to a previously built configfile, run:
```
./run_elvers config.yaml preprocess --print_params
```

A small set of parameters should print to your console:

```
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

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md)
