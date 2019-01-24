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

This is the default workflow for `eelpond`. To run:
```
./run_eelpond nema-test.yaml
# or
./run_eelpond nema-test.yaml full
```
This will run a small set of Nematostella vectensis test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16)).

## Running Your Own Data

To run your own data, you'll need to create two files, a tsv file containing your sample info, and a yaml file containing basic configuration info. 

IMPORTANT: The sample info must be in a **properly formatted** tsv file. The easiest way to do this is to copy the test data tsv and modify:
```
cp nema_samles.tsv my-samples.tsv
```
Now modify  `my-samples.tsv` with your sample information.

Next, build a configfile to edit:
```
./run_eelpond my-config.yaml --build_config
```
This configfile will contain all the default parameters for each step of the pipeline you target. If you don't specify any targets, it will default to the "full" Eel Pond Protocol pipeline, which executes read preprocessing, assembly, annotation, and quantification.

Please see the documentation file for each individual program (linked above) for what parameters to modify.

The configfile should look something like this:
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

