# Eel Pond Protocol Workflow

The Eel Pond protocol (which inspired the eelpond name) included line-by-line commands that the user could follow along with using a test dataset provided in the instructions. We have re-implemented the protocol here to enable automated de novo transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command. See more about this protocol [here](https://eel-pond.readthedocs.io/en/latest/).

The "Eel Pond" Protocol for RNAseq consists of:

  - [trimmomatic](trimmomatic.md) adapter and read quality trimming
  - [fastqc](fastqc.md) read qc evaluation
  - [khmer](khmer.md) k-mer trimming and (optional) digital normalization
  - [trinity](trinity.md) *de novo* assembly
  - [dammit](dammit.md) annotation
  - [salmon](salmon.md) read quantification to the trinity assembly
  - [deseq2](deseq2.md) differential expression analysis

## Running Test Data

This is the default workflow. To run:
```
./run_eelpond nema-test.yaml
# or
#./run_eelpond nema-test.yaml default
```
This will run a small set of Nematostella vectensis test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16)).

## Running Your Own Data

Set sample info and build a configfile first (see [Understanding and Configuring Workflows](about_and_configure.md)).

To build a config, run:

```
./run_eelpond ep.yaml --build_config
```

The resulting `ep.yaml` configfile for this workflow will look something like this. The order of the parameters may be different and does not affect the order in which steps are run. Please see the documentation file for each individual program (linked above) for what parameters to modify.

```
  ####################  Eelpond Pipeline Configfile  ####################
basename: eelpond
experiment: _experiment1
samples: samples.tsv ### PATH TO YOUR SAMPLE FILE GOES HERE

  ####################  assemble  ####################
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
trimmomatic:
  adapter_file:
    pe_name: TruSeq3-PE.fa
    pe_url: https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa
    se_name: TruSeq3-SE.fa
    se_url: https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa
  extra: ''
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
trinity:
  add_single_to_paired: false
  extra: ''
  input_kmer_trimmed: true
  input_trimmomatic_trimmed: false
  max_memory: 30G
  seqtype: fq

  ####################  annotate  ####################
dammit:
  busco_group:
  - metazoa
  - eukaryota
  db_dir: databases
  db_extra: ''
sourmash:
  extra: ''

  ####################  quantify  ####################
salmon:
  input_trimmomatic_trimmed: True
  index_params:
    extra: ''
  quant_params:
    extra: ''
    libtype: A

  ####################  diffexp  ####################
deseq2:
  contrasts:
    time0-vs-time6:
    - time0
    - time6
  gene_trans_map: true
  pca:
    labels:
    - condition
```

