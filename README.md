# eelpond

[![Build Status](https://travis-ci.org/dib-lab/eelpond.svg?branch=master)](https://travis-ci.org/dib-lab/eelpond)


```
                           ___
                        .-'   `'.
                       /         \
                      |           ;
                      |           |           ___.--,
             _.._     |O)  ~  (O) |    _.---'`__.-( (_.       
      __.--'`_.. '.__.\      '--. \_.-' ,.--'`     `""`
     ( ,.--'`   ',__ /./;     ;, '.__.'`    __
     _`) )  .---.__.' / |     |\   \__..--""  """--.,_
    `---' .'.''-._.-'`_./    /\ '.  \_.-~~~````~~~-.__`-.__.'
          | |  .' _.-' |    |  \  \  '.
           \ \/ .'     \    \   '. '-._)
            \/ /        \    \    `=.__`-~-.
            / /\         `)   )     / / `"".`\
      , _.-'.'\ \        /   /     (  (   /  /
       `--~`  )  )    .-'  .'       '.'. |  (
             (/`     (   (`           ) ) `-;
              `       '--;            (' 

```

eelpond is a transcriptome assembly, annotation and differential
expression workflow system for Illumina RNAseq. The default workflow
includes Trinity, dammit, and deseq2.

eelpond started as a snakemake update of the Eel Pond Protocol for *de
novo* RNAseq analysis. This assembly and annotation protocol was used
for the
[re-assembly, quality evaluation, and annotation of 678 microbial eukaryotic reference transcriptomes](https://academic.oup.com/gigascience/advance-article/doi/10.1093/gigascience/giy158/5241890),
among others.

eelpond has evolved to include workflows for
(mostly) RNA data, which can all be run via the `eelpond` workflow
wrapper. `eelpond` uses [snakemake](https://snakemake.readthedocs.io)
for workflow management and [conda](https://conda.io/docs/) for
software installation. The code can be found
[here](https://github.com/dib-lab/eelpond).


## Getting Started

Linux is the recommended OS. Nearly everything also works on MacOSX, but some programs (fastqc, Trinity) are troublesome.

If you don't have conda yet, install [miniconda](https://conda.io/miniconda.html) (for Ubuntu 16.04 [Jetstream image](https://use.jetstream-cloud.org/application/images/107)):
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
echo export PATH="$HOME/miniconda3/bin:$PATH" >> ~/.bash_profile
source ~/.bash_profile
```

Now, get the eelpond code
```
git clone https://github.com/dib-lab/eelpond.git
cd eelpond
```

Create a conda environment with all the dependencies for eelpond
```
conda env create --file environment.yml -n eelpond
```

Activate that environment. You'll need to do this anytime you want to run eelpond
```
conda activate eelpond
```
Now you can start running workflows on test data!

## Default workflow: Eel Pond Protocol for *de novo* RNAseq analysis

The Eel Pond protocol (which inspired the `eelpond` name) included line-by-line commands that the user could follow along with using a test dataset provided in the instructions. We have re-implemented the protocol here to enable automated *de novo* transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command. See more about this protocol [here](eel_pond_workflow.md).

To run the default workflow (assembly and annotation), do:
```
./run_eelpond examples/nema.yaml default
```
This will download, assemble, and annotate a small set of _Nematostella vectensis_ test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16)).

Note: this will take about 30 minutes the first time you run it,
because snakemake needs to install the various software packages.

## Running Your Own Data

To run your own data, you'll need to create two files:

  - a `tsv` file containing your sample info
  - a `yaml` file containing basic configuration info. This file must specify read inputs (path to `tsv`) or assembly inputs (via `assemblyinput`)

Generate these by following instructions here: [Understanding and Configuring Workflows](https://dib-lab.github.io/eelpond/configure).


## Available Workflows

Currently, all workflows require a properly-formatted read inputs `tsv` file as input. Some workflows, e.g. `annotation` can work on either on a _de novo_ transcriptome or on previously-generated assemblies. To add an assembly as input, specify it via `assemblyinput` in the `yaml` config file, as described in [Understanding and Configuring Workflows](https://dib-lab.github.io/eelpond/configure). 


**workflows**

  - preprocess: Read Quality Trimming and Filtering (fastqc, trimmomatic)
  - kmer_trim: Kmer Trimming and/or Digital Normalization (khmer)
  - assemble: Transcriptome Assembly (trinity)
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
./run_eelpond examples/nema.yaml --print_workflows
```

Each included tool can also be run independently, if appropriate input files are provided. This is not always intuitive, so please see our documentation for running each tools for details (described as "Advanced Usage"). To see all available tools, run:

```
./run_eelpond examples/nema.yaml --print_rules
```

## Additional Info

See the help, here:
```
./run_eelpond -h
```

**References:**  

  * [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
  * [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
  * [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
  * [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)


