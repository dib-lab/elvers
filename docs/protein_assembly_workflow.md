# Protein Assembly Workflow

The protein assembly workflow relies on the PLASS assembler and some downstream protein mapping tools.

It consists of:  

  - [fastqc](fastqc.md)
  - [trimmomatic](trimmomatic.md)
  - [khmer](khmer.md)
  - [PLASS](plass.md)
  - [pear](pear.md)
  - [paladin](paladin.md)

## Running Test Data

```
./run_eelpond examples/nema.yaml protein_assembly
```
This will run a small set of Nematostella vectensis test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16)).

## Running Your Own Data

Set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To build a config, run:

```
./run_eelpond prot.yaml protein_assembly --build_config
```

The resulting `prot.yaml` configfile for this workflow will look something like this. The order of the parameters may be different and does not affect the order in which steps are run. Please see the documentation file for each individual program (linked above) for what parameters to modify.

```
  ####################  Eelpond Pipeline Configfile  ####################
basename: eelpond
experiment: _experiment1
samples: samples.tsv

  ####################  protein_assembly  ####################
fastqc:
  extra: ''
get_data:
  download_data: false
khmer:
  C: 3
  Z: 18
  coverage: 20
  diginorm: true
  extra: ''
  ksize: 20
  memory: 4e9
paladin:
  alignment_params:
    extra: ''
    f: 125
  index_params:
    reference_type: '3'
pear:
  extra: ''
  input_kmer_trimmed: false
  input_trimmomatic_trimmed: true
  max_memory: 4G
  pval: 0.01
plass:
  add_single_to_paired: false
  extra: ''
  input_kmer_trimmed: true
  input_trimmomatic_trimmed: false
sourmash:
  extra: ''
  k_size: 31
  scaled: 1000
trimmomatic:
  adapter_file:
    pe_path: ep_utils/TruSeq3-PE-2.fa
    se_path: ep_utils/TruSeq3-SE.fa
  extra: ''
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25
```
