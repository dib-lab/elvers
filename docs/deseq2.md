# Differential expression analysis with DESeq2

We can use [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to compare gene expression differences in samples between experimental conditions.

## Quickstart: Running DESeq2 via elvers

We recommend you run `deseq2` via the [diffexp](diffexp.md) subworkflow.

If you want to run it as a standalone program instead, you need to have generated read quantification data via `salmon`.

   1) If you have salmon results, run:
   
   ```
   elvers examples/nema.yaml deseq2
   ```

   2) If not, you need to run `salmon` and any other missing steps. It's probably best to run the [diffexp](diffexp.md) subworkflow, but you can also try: 

   ```
   elvers examples/nema.yaml salmon deseq2
   ```

## DESeq2 Commands

This pipeline uses snakemake to run a few R scripts to conduct basic differential expression analysis. We read in transcript abundance information (generated with [salmon](salmon.md)) via [tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html). Note that in the salmon step, we combine files of all "units" within a sample in order to then conduct differential expression at the sample level. 

We assume the assembly has a gene-to-transcript map, such as the one produced via trinity. This is a tab separated file (`transcript \t gene`) that enables count data to be aggregated at the gene level prior to differntial expression analysis. This is recommended, see [Soneson et al, 2016](https://f1000research.com/articles/4-1521/v2). However, if you do not have this mapping, we provide an option to conduct differential expression at the transcript level via the config (see "Customizing DESeq2 Parameters" section, below).

After reading in count data, we take in two additional pieces of information: first, the sample names in the `samples.tsv` document, and second the desired `contrast`, provided as part of the DESeq2 parameters, below. We store all data in an `.rds` r data format to support easy reloading of this data for additional user analyses. In addition, we plot a PCA of the normalized counts and perform a standard DESeq2 analysis and print a `tsv` of results for each contrast specified in
the deseq2 params.

You can find these R scripts in the `elvers` [github repo](https://github.com/dib-lab/elvers/tree/master/rules/deseq2). The snakemake rules and scripts were modified from [rna-seq-star-deseq2 workflow](https://github.com/snakemake-workflows/rna-seq-star-deseq2) and our own data analysis and workshops, e.g.[DIBSI-RNAseq](https://dibsi-rnaseq.readthedocs.io/en/latest/DE.html). 

## Modifying Params for DESeq2

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `deseq2` rule, run
```
elvers config deseq2 --print_params
```
This will print the following:
```
  ####################  deseq2  ####################
deseq2:
  contrasts:
    time0-vs-time6:
    - time0
    - time6
  gene_trans_map: true
  pca:
    labels:
    - condition
  #####################################################
```
The default `contrasts` reflect the `condition` information in the test data `nema_samples.tsv`. Please modify the contrasts to the reflect your data. Multiple contrasts should be supported: each contrast needs a name, and a list below it specifying the conditions to compare, e.g.:
```
  contrasts:
    my-contrast:
      - conditionA
      - conditionB
```
The `pca labels` should not be changed unless you need to change the name of the `condition` column in the `samples.tsv`. This functionality hasn't been extensively tested, so file an issue if something goes wrong!

Be sure the modified lines go into the config file you're using to run `elvers` (see [Understanding and Configuring Workflows](configure.md)).


## References

  * [Documentation for DESeq2 with example analysis](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
  * [Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
  * [Love et al. 2016](https://www.nature.com/nbt/journal/v34/n12/full/nbt.3682.html)

Additional links:

  * [DE lecture by Jane Khudyakov, July 2017](https://github.com/Open-Data-Science-at-SIO/RNAseq-workshop-2017/blob/master/_static/Jane_differential_expression.pdf)
  * [Example DE analysis from two populations of killifish! (Fundulus heteroclitus MDPL vs. MDPL)](http://htmlpreview.github.io/?https://github.com/ljcohen/Fhet_MDPL_MDPP_salinity_DE/blob/master/Fhet_MDPL_v_MDPP_interactiononly_FW_BW.html)
  * [A Review of Differential Gene Expression Software for mRNA sequencing](https://github.com/ljcohen/ECE221_final_project/blob/master/Cohen_Li_ECE221_review-differential-gene.pdf)


## Snakemake Rules

For snakemake afficionados, see the deseq2 rules on [github](https://github.com/dib-lab/elvers/blob/master/rules/deseq2/deseq2.rule).



