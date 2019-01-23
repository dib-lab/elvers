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
Snakemake update of the Eel Pond Protocol for *de novo* RNAseq analysis


**OSX issues:**
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
conda env create --file environment.yml -n eelpond
```

Activate that environment. You'll need to do this anytime you want to run eelpond
```
conda activate eelpond
```

Now, grab the test data, and untar:
```
curl -L https://osf.io/chb7z/download -o nema_testdata.tar.gz
tar xvf nema_testdata.tar.gz
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



