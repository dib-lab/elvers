# EelPond

# Snakemake update of the Eel Pond Protocol for *de novo* RNAseq analysis

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
This is a lightweight protocol for assembling up to a few hundred million mRNAseq reads, annotating the resulting assembly, and doing differential expression analysis. The input is short-insert paired-end Illumina reads. This protocol can be run in a single command because it uses the snakemake automated workflow management system.

Previous versions of this protocol included line-by-line commands that the user could follow along with using a test dataset provided in the instructions. Since the recent development of [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management tool and [snakemake-wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/) to manage sofware installation of commonly-used bioinformatics tools, we have re-implemented the Eel Pond Protocol to make it easier for users to install software and run a *de novo* transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command.

The software for this protocol can be found [here](https://github.com/dib-lab/eelpond). 


# Getting Started

At the moment, only Linux is supported. OSX issues:
  - fastqc fails about half the time
  - Trinity assembler does not work


Install [miniconda](https://conda.io/miniconda.html) (for Ubuntu 16.04 [Jetstream image](https://use.jetstream-cloud.org/application/images/107)):
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
conda env create --file ep_utils/eelpond_environment.yaml -n eelpond
```

Activate that environment. You'll need to do this anytime you want to run eelpond
```
source activate eelpond
```

Now you can start running eelpond!

To test a "full" workflow, consisting of read pre-processing, kmer trimming, Trinity assembly, dammit annotation and salmon quantification:
```
./run_eelpond nema-test full
```
These will run a small set of _Nematostella vectensis_ test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16))


You can also run individual tools or subworkflows independently:
```
./run_eelpond nema-test preprocess
./run_eelpond nema-test trimmomatic
```

See the help, here:
```
./run_eelpond -h
```

**Running your own data:**

To run your own data, you'll need to create two files, a `tsv` file containing 
your sample info, and a `yaml` file containing basic configuration info. To start,
copy the test data files so you can modify them.

```
cp nema_samples.tsv <my-tsv-name.tsv>
```

Next, build a configfile to edit:

```
./run_eelpond config_name --build_config

```
This configfile will contain all the default paramters for each step of the pipeline you target.
If you don't specify any targets, it will default to the "full" pipeline, which executes read
preprocessing, assembly, annotation, and quantification.

Then, modify this configfile as necessary. 
The essential component is the `samples.tsv` file, which points `eelpond` to your sample files.


**References:**  

  * [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
  * [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
  * [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
  * [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)


**available workflows:**  

  - preprocess: Read Quality Trimming and Filtering (fastqc, trimmomatic)
  - kmer_trim: Kmer Trimming and/or Digital Normalization (khmer)
  - assemble: Transcriptome Assembly (trinity)
  - assemblyinput: Specify assembly for downstream steps
  - annotate : Annotate the transcriptome (dammit, sourmash)
  - quantify: Quantify transcripts (salmon) 
  - full: preprocess, kmer_trim, assemble, annotate, quantify 



