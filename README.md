# 2018-snakemake-eel-pond
Snakemake update of the Eel Pond Protocol for *de novo* RNAseq analysis

Install [miniconda](https://conda.io/miniconda.html) :

For Ubuntu 16.04 ([Jetstream image](https://use.jetstream-cloud.org/application/images/107))
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
git clone https://github.com/dib-lab/2018-snakemake-eel-pond.git
cd 2018-snakemake-eel-pond
snakemake --use-conda --configfile cfp.yml
```


**References:** 
* [eel-pond protocol docs](http://eel-pond.readthedocs.io/en/latest/)
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


*snakemake style follows [rna-seq-star example workflow](https://github.com/snakemake-workflows/rna-seq-star-deseq2)*
