# Eel Pond

This is a lightweight protocol for assembling up to a few hundred million mRNAseq reads, annotating the resulting assembly, and doing differential expression analysis. The input is short-insert paired-end Illumina reads. This protocol can be run in a single command because it uses the snakemake automated workflow management system.

Previous versions of this protocol included line-by-line commands that the user could follow along with using a test dataset provided in the instructions. Since the recent development of [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management tool and [snakemake-wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/) to manage sofware installation of commonly-used bioinformatics tools, we have re-implemented the Eel Pond Protocol to make it easier for users to install software and run a *de novo* transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command.

The software for this protocol can be found [here](https://github.com/dib-lab/eelpond). 

To run the protocol on your own computer system (requires Ubuntu 16.04):

Install [miniconda](https://conda.io/miniconda.html) (for Ubuntu 16.04 [Jetstream image](https://use.jetstream-cloud.org/application/images/107)).
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
echo export PATH="$HOME/miniconda3/bin:$PATH" >> ~/.bash_profile
source ~/.bash_profile
```

Then,

```
conda install -c bioconda -c conda-forge -y snakemake
```

And:

```
git clone https://github.com/dib-lab/eelpond.git
cd eelpond
snakemake --use-conda --configfile cfp.yml
```

**References:** 
* [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
* [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
* [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
* [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)

**intended workflows:**
  - [Read Quality Trimming and Filtering](QC.md)
  - [Digital Normalization](Diginorm.md)
  - [Assembly](Assembly.md)
  - [Quality Assessment](Quality.md)
  - [Annotation](Annotation.md)
  - [Transcript Quantification](Quant.md)
  - [Differential Expression](DE.md)


*snakemake style follows [rna-seq-star example workflow](https://github.com/snakemake-workflows/rna-seq-star-deseq2)*
