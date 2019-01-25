# Mapping to Protein Assemblies using Paladin

We use `Paladin` to map back to protein assemblies ([plass](plass.md)). `Paladin` solves the problem of mapping *nucleotide* sequences back to *amino acid* sequences.

From the [documentation](https://github.com/twestbrookunh/paladin):

Protein ALignment And Detection INterface

PALADIN is a protein sequence alignment tool designed for the accurate functional characterization of metagenomes.

PALADIN is based on BWA, and aligns sequences via read-mapping using BWT. PALADIN, however, offers the novel approach of aligning in the protein space. During the index phase, it processes the reference genome's nucleotide sequences and GTF/GFF annotation containing CDS entries, first converting these transcripts into the corresponding protein sequences, then creating the BWT and suffix array from these proteins. The process of translatation is skiped when providing a protein reference file (e.g., UniProt) for mapping. During the alignment phase, it attempts to find ORFs in the read sequences, then converts these to protein sequences, and aligns to the reference protein sequences.

PALADIN currently only supports single-end reads (or reads merged with FLASH, PEAR, abyss-mergepairs), and BWA-MEM based alignment. It makes use of many BWA parameters and is therefore compatible with many of its command line arguments.

PALADIN may output a standard SAM file, or a text file containing a UniProt-generated functional profile. This text file may be used for all downstream characterizations.

## Quickstart: Running Paladin with eelpond

We recommend you run `paladin` as part of the [Protein Assembly](protein_assembly_workflow.md) or [paladin_map](paladin_map.md) workflows.

```
./run_eelpond nema-test protein_assembly
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

## Modifying Params for Paladin:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](about_and_configure.md)).

To see the available parameters for the `Paladin` rule, run
```
./run_eelpond config paladin --print_params
```
This will print the following:
```
  ####################  paladin  ####################
paladin:
  alignment_params:
    extra: ''
    f: 125
  index_params:
    reference_type: '3'
    gff_file: '' 
  #####################################################
```
In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass in additional parameters to `paladin align`,  e.g.:
```
  extra: ' --some_param that_param '
```
Please see the [Paladin documentation](https://github.com/twestbrookunh/paladin) for info on the params.

Be sure the modified lines go into the config file you're using to run `eelpond` (see [Understanding and Configuring Workflows](about_and_configure.md)).

## Advanced Usage: Running PALADIN as a standalone rule

You can run paladin as a standalone rule, instead of withing a larger `eelpond` workflow. However, to do this, you need to make sure the input files are available.

For paladin, you need both 1) an assembly, and 2) trimmed (and merged) input files. The assembly can be generated via another workflow, or passed to `eelpond` via the configfile.

Specifying an assembly:
    1) If you've alread run read trimming and want to use a Trinity assembly generated via `eelpond`, you can run:
    ```
    ./run_eelpond my_config plass_assemble paladin
    ```
    If you've already run the assembly, `eelpond` will just use this info to locate that assembly.

    2) Alternatively, you can input an assembly via the [assemblyinput](assemblyinput.md) utility rule:
    ```
    ./run_eelpond assemblyinput paladin
     ```
    with an assembly in your `yaml` configfile, e.g.:
    ```
    assemblyinput:
      assembly: rna_testdata/nema.fasta
      gene_trans_map:  rna_testdata/nema.fasta.gene_trans_map #optional
      assembly_extension: '_plass'
    ```
    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have a gene to transcript map, please specify it as well. If not, delete this line from    your `config`. The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Make sure you set the `assembly_extension` parameter to plass, as paladin only works on plassassemblies (for the
    moment). **Note: Please don't use additional underscores (`_`) in this extension!**. For more details, see the [assemblyinput documentation](assemblyinput.md).

Specifying input reads:

If you haven't yet run read trimming and merging, you'll also need to run those steps:
```
./run_eelpond my_config get_data trimmomatic pear paladin
```
with one of the options to specify an assembly (above).


## PALADIN eelpond rule

We wrote new snakemake wrappers for [paladin index](https://github.com/dib-lab/eelpond/blob/master/rules/paladin/paladin-index.py) and [paladin align](https://github.com/dib-lab/eelpond/blob/master/rules/paladin/paladin-align.py) to run PALADIN via snakemake. These wrappers have not been added to the official repo yet, but feel free to use as needed.

For snakemake afficionados, see our paladin rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/paladin/paladin.rule).

## Citation

If you use PALADIN, please cite [Westbrook et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5423455/).
