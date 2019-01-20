# Mapping to Protein Assemblies using Paladin

We use `Paladin` to map back to protein assemblies ([plass](plass.md)). `Paladin` solves the problem of mapping *nucleotide* sequences back to *amino acid* sequences.

From the [documentation](https://github.com/twestbrookunh/paladin):

Protein ALignment And Detection INterface

PALADIN is a protein sequence alignment tool designed for the accurate functional characterization of metagenomes.

PALADIN is based on BWA, and aligns sequences via read-mapping using BWT. PALADIN, however, offers the novel approach of aligning in the protein space. During the index phase, it processes the reference genome's nucleotide sequences and GTF/GFF annotation containing CDS entries, first converting these transcripts into the corresponding protein sequences, then creating the BWT and suffix array from these proteins. The process of translatation is skiped when providing a protein reference file (e.g., UniProt) for mapping. During the alignment phase, it attempts to find ORFs in the read sequences, then converts these to protein sequences, and aligns to the reference protein sequences.

PALADIN currently only supports single-end reads (or reads merged with FLASH, PEAR, abyss-mergepairs), and BWA-MEM based alignment. It makes use of many BWA parameters and is therefore compatible with many of its command line arguments.

PALADIN may output a standard SAM file, or a text file containing a UniProt-generated functional profile. This text file may be used for all downstream characterizations.

## Quickstart: Running Paladin with eelpond

We recommend you run `paladin` as part of the `paladin_map` pipeline

```
./run_eelpond nema-test paladin_map
```
This will run trimmomatic trimming prior to PEAR merging (for paired end reads) and paladin mapping. 

Note, `PEAR` only works on paired end reads, as there's nothing to do for single end reads. By default, SE reads are taken directly from `trimmomatic` output.

## Paladin Commands

Eelpond first indexes the protein assembly, and then maps reads back. On the command line, these commands would look like this: 
```
paladin index 3 index_basename \ 
  gff_annotation extra
```
We have not yet added dammit annotation gff at this step, but this could be done :). You can also add any `gff3` annotation via the params.

```
paladin align -f 250 -t snakemake.threads \
  index_basename reads_file | samtools view -Sb - > output.bam
```

## Customizing Paladin parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a Paladin configfile you can modify, run:
```
./run_eelpond paladin.yaml paladin --build_config
```
The output should be a small `yaml` configfile that contains:
```
  ####################  paladin  ####################
paladin:
  alignment_params:
    extra: ''
    f: 125
  index_params:
    reference_type: '3'
    gff_file: ''
```
In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass in additional parameters to `paladin align`,  e.g.:

```
  extra: ' --some_param that_param '
```
Please see the [Paladin documentation](https://github.com/twestbrookunh/paladin) for info on the params.

Override default params by modifying any of these lines, and placing them in the config file you're using to run `eelpond`. Here, we just generated params for `PALADIN`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file. 

For example, run the following to build default parameters for `trimmomatic`, `fastqc`, `khmer`, `plass`, `pear`, and `paladin`, which you can then edit in `my-workflow.yaml`:
```
./run_eelpond my-workflow.yaml plass_assemble paladin_map --build_config
``` 

## PALADIN eelpond rule

We wrote new snakemake wrappers for [paladin index](https://github.com/dib-lab/eelpond/blob/master/rules/paladin/paladin-index.py) and [paladin align](https://github.com/dib-lab/eelpond/blob/master/rules/paladin/paladin-align.py) to run PALADIN via snakemake. These wrappers have not been added to the official repo yet, but feel free to use as needed.

For snakemake afficionados, see our paladin rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/paladin/paladin.rule).

## Citation

If you use PALADIN, please cite [Westbrook et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5423455/).
