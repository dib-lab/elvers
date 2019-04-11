# elvers

[![Build Status](https://travis-ci.org/dib-lab/elvers.svg?branch=master)](https://travis-ci.org/dib-lab/elvers)


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
elvers started as a snakemake update of the Eel Pond Protocol for *de novo* RNAseq analysis. It has evolved slightly to enable a number of workflows for (mostly) RNA data, which can all be run via the `elvers` workflow wrapper. `elvers` uses [snakemake](https://snakemake.readthedocs.io) for workflow management and [conda](https://conda.io/docs/) for software installation. The code can be found [here](https://github.com/dib-lab/elvers).


## Getting Started

Linux is the recommended OS. Nearly everything also works on MacOSX, but some programs (fastqc, Trinity) are troublesome.

If you don't have conda yet, install [miniconda](https://conda.io/miniconda.html) (for Ubuntu 16.04 [Jetstream image](https://use.jetstream-cloud.org/application/images/107)):
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Be sure to answer 'yes' to all yes/no questions. You'll need to restart your terminal for conda to be active.

## Create a working environment and install `elvers`!

`elvers` needs a few programs installed in order to run properly. To handle this, we run `elvers` within a conda environment that contains all dependencies.

Get the **`elvers`** code
```
git clone https://github.com/dib-lab/elvers.git
cd elvers
```

When you first get **`elvers`**, you'll need to create this environment on your machine:
```
conda env create --file environment.yml -n elvers-env
```

Now, activate that environment:
```
conda activate elvers-env
```
To deactivate after you've finished running `elvers`, type `conda deactivate`. You'll need to reactivate this environment anytime you want to run elvers.

Now. install the `elvers` package. 
```
pip install -e .
```
Now you can start running workflows on test data!

## Default workflow: Eel Pond Protocol for *de novo* RNAseq analysis

The Eel Pond protocol (which inspired the `elvers` name) included line-by-line commands that the user could follow along with using a test dataset provided in the instructions. We have re-implemented the protocol here to enable automated *de novo* transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command. See more about this protocol [here](eel_pond_workflow.md).

To test the default workflow:
```
elvers examples/nema.yaml default
```
This will download and run a small set of _Nematostella vectensis_ test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16))

## Running Your Own Data

To run your own data, you'll need to create one or more files:

  - a `yaml` file containing basic configuration info
 
 This `yaml` config file must specify **either**:
   - a `tsv` file containing your read sample info
   - a reference file input (`.fasta` file and optional `gene_trans_map`)


Generate these files by following instructions here: [Understanding and Configuring Workflows](https://dib-lab.github.io/elvers/configure).


## Available Workflows

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
elvers examples/nema.yaml --print_workflows
```

Each included tool can also be run independently, if appropriate input files are provided. This is not always intuitive, so please see our documentation for running each tools for details (described as "Advanced Usage"). To see all available tools, run:

```
elvers examples/nema.yaml --print_rules
```

## Citation information
  
This is *pre-publication* code; a manuscript is in preparation. Please contact the authors for the current citation information if you wish to use it and cite it.


## Additional Info

See the help, here:
```
elvers -h
```

**References:**  

  * [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
  * [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
  * [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
  * [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)


