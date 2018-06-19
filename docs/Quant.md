# Quantification with Salmon

We will use [Salmon](http://salmon.readthedocs.org/en/latest/) to
quantify expression. Salmon is a new breed of software for quantifying RNAseq reads that is both really fast and takes
transcript length into consideration ([Patro et al. 2015](http://dx.doi.org/10.1038/nmeth.4197)).

For further reading, see

  * Intro blog post: http://robpatro.com/blog/?p=248
  * A 2016 blog post evaluating and comparing methods [here](https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/)
  * Salmon github repo [here](https://github.com/COMBINE-lab/salmon)
  * https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/salmon.rst
  * http://angus.readthedocs.io/en/2016/rob_quant/tut.html
  * https://2016-aug-nonmodel-rnaseq.readthedocs.io/en/latest/quantification.html

The two most interesting files are `salmon_quant.log` and
`quant.sf`. The latter contains the counts; the former contains the
log information from running things.

We recommend quantifying using the Trinity transcriptome assembly fasta file, which will give expression values for each contig, like this in `quant.sf`:
```
Name                  Length    EffectiveLength    TPM    NumReads
TRINITY_DN2202_c0_g1_i1    210    39.818    2.683835    2.000000
TRINITY_DN2270_c0_g1_i1    213    41.064    0.000000    0.000000
TRINITY_DN2201_c0_g1_i1    266    69.681    0.766816    1.000000
TRINITY_DN2222_c0_g1_i1    243    55.794    2.873014    3.000000
TRINITY_DN2291_c0_g1_i1    245    56.916    0.000000    0.000000
TRINITY_DN2269_c0_g1_i1    294    89.251    0.000000    0.000000
TRINITY_DN2269_c1_g1_i1    246    57.479    0.000000    0.000000
TRINITY_DN2279_c0_g1_i1    426    207.443    0.000000    0.000000
TRINITY_DN2262_c0_g1_i1    500    280.803    0.190459    1.000912
TRINITY_DN2253_c0_g1_i1    1523    1303.116    0.164015    4.000000
TRINITY_DN2287_c0_g1_i1    467    247.962    0.000000    0.000000
TRINITY_DN2287_c1_g1_i1    325    113.826    0.469425    1.000000
TRINITY_DN2237_c0_g1_i1    306    98.441    0.542788    1.000000
TRINITY_DN2237_c0_g2_i1    307    99.229    0.000000    0.000000
TRINITY_DN2250_c0_g1_i1    368    151.832    0.000000    0.000000
TRINITY_DN2250_c1_g1_i1    271    72.988    0.000000    0.000000
TRINITY_DN2208_c0_g1_i1    379    162.080    1.978014    6.000000
TRINITY_DN2277_c0_g1_i1    269    71.657    0.745677    1.000000
TRINITY_DN2231_c0_g1_i1    209    39.409    0.000000    0.000000
TRINITY_DN2231_c1_g1_i1    334    121.411    0.000000    0.000000
TRINITY_DN2204_c0_g1_i1    287    84.121    0.000000    0.000000
```
