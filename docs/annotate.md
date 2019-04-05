# Annotate Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting withing `elvers`. The "annotate" subworkflow annotates a transcriptome. It requires an assembly to be provided, either by running an assembly or providing one in your configfile. 

At the moment, this workflow consists of:
 
  - [dammit](dammit.md)

## Quickstart

If you're starting a new `elvers` run using a reference file (even if you have previously built a _de novo_ assembly via `elvers`), you need to help `elvers` find that assembly.


    Scenario 1: You're starting from your own reference file:

    You need to provide the reference in your `config.yaml` file:
   
    ```
    get_reference:
      reference: input reference fasta REQUIRED
      gene_trans_map: OPTIONAL: provide a gene to transcript map for input transcriptome
      reference_extension: '_input' OPTIONAL, changes naming
      download_ref: the reference entry above is a link that needs to be downloaded
      use_ftp: download via ftp instead of http
    ```
    Once you add this to your configfile (e.g. `my_config.yaml`), you can run a reference/assembly-based workflow such as annotate.
    ```
    elvers my_config.yaml annotate
    ```
    For more details on reference specification, see the [get_reference documentation](get_reference.md). For annotation configuration, see below.

    Scenario 2: You've previously run an assembly program via `elvers`:

    You _can_ provide the built reference in the same manner as above. However, if you're running more workflows in the same directory, you can also just specify the name of the assembly program that you used to generate the assembly. This will *not* rerun the assembly (unless you provide new input files). Instead, this will allow `elvers` to know where to look for your previously-generated reference file. Because we enable multiple referene generation programs, we don't want to assume
   which reference you'd like to use for downstream steps (in fact, if you provide multiple references, `elvers` will run the downstream steps on all references, assuming you provide unique `reference_extension` parameters so that the references are uniquely named.
   
    Example: You've previously run the trinity assembly, and want to annotate it.
   
    ```
    elvers examples/nema.yaml assemble annotate
    ```
    Here, the `assemble` workflow just enables `elvers` to locate your assembly file for `annotate`.  

   
## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Dammit will output files in the `annotation` subdirectory of this output directory. The annotated fasta file will be `ASSEMBLY.dammit.fasta` and the annotation `gff3` file will be `ASSEMBLY.dammit.gff3`. Dammit will also produce a number of intermediate files that will be contained within an `ASSEMBLY.fasta.dammit` folder.

## Configuring the annotate subworkflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](configure.md).

If you want to add the `annotate` program parameters to a previously built configfile, run:
```
elvers config.yaml annotate --print_params
```

A small set of parameters should print to your console:

```
 ####################  annotate  ####################
dammit:
  busco_group:     # specify all busco groups below here
  - metazoa
  - eukaryota
  db_dir: databases   # specify location for databases (or previously-installed databases)
  db_install_only: False   # just install databases, don't run annotation
  db_extra:
  annot_extra: ' --quick '
  #######################################################
```

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [dammit](dammit.md)
