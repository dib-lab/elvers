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
conda install -c bioconda -c conda-forge -y snakemake
```


Run:

```
git clone https://github.com/dib-lab/eelpond.git
cd eelpond
snakemake --use-conda --configfile rna_testdata/nema_config.yaml
```


**References:**  

  * [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
  * [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
  * [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
  * [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)


**intended workflows:**  

  - Read Quality Trimming and Filtering
  - Digital Normalization
  - Assembly
  - Quality Assessment
  - Annotation
  - Transcript Quantification 
  - Differential Expression


*snakemake style builds off of [rna-seq-star example workflow](https://github.com/snakemake-workflows/rna-seq-star-deseq2)*
