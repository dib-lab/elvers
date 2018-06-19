# Evaluating your transcriptome assembly

We will be using Transrate and BUSCO!

## BUSCO

* **B**enchmarking **U**niversal **S**ingle **C**opy **O**rthologs
* Eukaryota database has 303 genes
* Metazoa database has 978 genes
* "Complete" lengths are within two standard deviations of the BUSCO group mean length
* Genes that make up the BUSCO sets for each major lineage are selected from orthologous groups with genes present as single-copy orthologs in at least 90% of the species. 

* Useful links:
  * Website with additional busco databases: [http://busco.ezlab.org/](http://busco.ezlab.org/)
  * Paper: [Simao et al. 2015](http://bioinformatics.oxfordjournals.org/content/31/19/3210)
  * [User Guide](http://gitlab.com/ezlab/busco/raw/master/BUSCO_v2.0_userguide.pdf)

Command:

```
run_BUSCO.py \
-i Trinity.fixed.fasta \
-o nema_busco_metazoa -l ~/busco/metazoa_odb9 \
-m transcriptome --cpu 2
```


## Transrate

[Transrate](http://hibberdlab.com/transrate/getting_started.html) serves two main purposes. It can compare two assemblies to see how similar they are. Or, it can give you a score which represents proportion of input reads that provide positive support for the assembly. We will use transrate to get a score for the assembly. Use the trimmed reads. For a further explanation of metrics and how to run the reference-based transrate, see the [documentation](http://hibberdlab.com/transrate/metrics.html) and the paper by [Smith-Unna et al. 2016](http://genome.cshlp.org/content/early/2016/06/01/gr.196469.115). 

* How do two transcriptomes compare with each other?

```
transrate --reference=Trinity.fixed.fasta --assembly=trinity-nematostella-raw.fa --output=full_v_subset
transrate --reference=trinity-nematostella-raw.fa --assembly=Trinity.fixed.fasta --output=subset_v_full
```
