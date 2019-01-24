# Kmer_trim Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting withing `eelpond`. The "kmer_trim" subworkflow conducts read quality trimming and kmer trimming. At the moment, this workflow consists of:
 
  - [get_data](get_data.md) - an `eelpond` utility
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md), run on both pre-trim and post-trim data
  - [khmer](khmer.md)


# Quickstart

To run the kmer_trim subworkflow, run: 

```
./run_eelpond nema-test kmer_trim
```

## Configuring the kmer_trim subworkflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](about_and_configure.md).

If you want to add the `kmer_trim` program parameters to a previously built configfile, run:
```
./run_eelpond config.yaml kmer_trim --print_params
```

A small set of parameters should print to your console:

```
 ####################  kmer_trim  ####################
get_data:
  download_data: false
  use_ftp: false
trimmomatic:
  adapter_file:
    pe_path: ep_utils/TruSeq3-PE-2.fa
    se_path: ep_utils/TruSeq3-SE.fa
  extra: ''
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25
fastqc:
  extra: ''
khmer:
  C: 3
  Z: 18
  coverage: 20
  diginorm: true
  extra: ''
  ksize: 20
  memory: 4e9
  #######################################################
```

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](about_and_configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [fastqc](fastqc.md)
  - [khmer](khmer.md)
