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

git submodule update --init --recursive #download test data submodule

#run eelpond
./run_eelpond nema-test
```

This will run a small set of _Nematostella vectensis_ test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16))

**Running your own data:**

To run your own data, you'll need to create two files, a `tsv` file containing 
your sample info, and a `yaml` file containing basic configuration info. To start,
copy the test data files so you can modify them.

```
cp nema-test.yaml <my_data_name.yaml>
cp nema_samples_v2.tsv <my-tsv-name.tsv>
```

Then, modify as necessary. Be sure to specify the new name of your `tsv` file within the `yaml` configuration file.


**References:**  

  * [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
  * [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
  * [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
  * [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)


**available workflows:**  

  - preprocess: Read Quality Trimming and Filtering (fastqc, trimmomatic)
  - ktrim: Kmer Trimming and/or Digital Normalization (khmer)
  - assembly: Transcriptome Assembly (trinity)
  - assembly_quality: Assess Assembly Quality (busco, sourmash, transrate) 
  - annotation : Annotate the transcriptome (dammit)
  - quantification: Quantify transcripts (salmon) 

  _in progress_
  - diffexp: Basic differential expression analyses (deseq2)


