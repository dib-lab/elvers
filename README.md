# 2018-snakemake-eel-pond
Snakemake update of the Eel Pond Protocol for *de novo* RNAseq analysis

run:
```
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


* migrated from [dahak-eel-pond](https://github.com/bluegenes/dahak-eel-pond) 
using ideas from [rna-seq-star example workflow](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
