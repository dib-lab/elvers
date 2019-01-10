# eelpond

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


* OSX issues:*
  - fastqc fails about half the time
  - Trinity assembler does not work


Install [miniconda](https://conda.io/miniconda.html) (for Ubuntu 16.04 [Jetstream image](https://use.jetstream-cloud.org/application/images/107)):
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
echo export PATH="$HOME/miniconda3/bin:$PATH" >> ~/.bash_profile
source ~/.bash_profile
```

```
conda install -c bioconda -c conda-forge -y snakemake yaml
```


Run Eelpond:

```
#get eelpond code
git clone https://github.com/dib-lab/eelpond.git
cd eelpond

# see the help:
./run_eelpond -h
```

If you're on Linux, run a full test:
```
#run test data
./run_eelpond nema-test full
```
This will run a small set of _Nematostella vectensis_ test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16))

If you're on OSX, Trinity will not work (for now). You can test all other steps:
individually:
```
./run_eelpond nema-test preprocess
```
or together:
```
./run_eelpond nema-test kmer_trim quantify annotate
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



