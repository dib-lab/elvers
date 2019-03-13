# Under the Hood

`eelpond` runs snakemake rules in workflows that are highly customizeable. The user chooses a tool or workflow, and `./run_eelpond` aggregates all of the default parameters for that tool or workflow, build end-stage target files, and passes this information into the snakefile. Snakemake then builds a DAG of the workflow and runs the workflow in an automated manner.

## Program Rules

The unit at the heart of all workflows is a standalone snakemake rule for a program function. Rules are all placed within subdirectories of the `rules` directory. For example, the `trinity` rule is found at `rules/trinity/trinity.rule`. Utility rules are instead found at, for example: `rules/utils/get_data.rule`. 

Each rule needs several files:

  - `rulename.rule`: the snakemake rule(s) for this program.
  - `rulename-env.yaml`: a conda environment file for this program. This is passed into the `conda` snakemake directive in the rule.
  - (optional) `rulename-wrapper.py` or `rulename-script.R`: a script that is used to more easily and reproducibly run the program.
  - `rulename_params.yaml`: a parameter file to direct program output and set default program parameters

The last file (`_params`) provides program defaults that can be overridden by user input from the main config file (via a nested dictionary update). The parameters do need to be exact, however. The `--print_params` and `--build_config` options in `./run_eelpond` are designed to facilitate this. 

## Rule Params

The rule params files (above) have two components: `eelpond_params` and `program_params`. `program_params` are exposed to the user when displaying parameters or building configfiles. These are relatively safe to modify, and most get passed into the program rule itself. On the other hand, `eelpond_params` are used internally to build the names of output files when they need to be passed to the Snakefile as "targets". These are never exposed to the user and should be set just once when a rule is built. 

For example, let's look at the `trinity.yaml` file:

```
trinity:
  eelpond_params:
    outdir: assembly
    extensions:
      assembly_extensions: # use this extension only for all output of this assembler
        - _trinity
      base:
        - .fasta
        - .fasta.gene_trans_map
  program_params:
    # input kmer-trimmed reads
    input_kmer_trimmed: True
    # input trimmed-reads
    input_trimmomatic_trimmed: False
    # do we want to assemble the single reads with pe reads?
    add_single_to_paired: False
    max_memory : 30G
    seqtype: fq
    extra: ''
```
The `program_params` allow the user to choose quality-trimmed or kmer-trimmed reads as input to trinity, and pass paramters for memory and sequence type. Whereever possible, program params also include an `extra` parameter that enable the user to pass any number of parameters directly to the command line of that program.

The `eelpond_params` specify that trinity produces two files, with the assembly_extension '_trinity'. These will be `BASENAME_trinity.fasta` and `BASENAME_trinity.fasta.gene_trans_map`. The `assembly_extension` allows us to have downstream steps that work on different assemblies, but trinity is an assembler, and should only ever produce an assembly with the extension `_trinity`. Assembly extensions are used downstream as well. For example, `paladin` only works on protein assemblies, and
thus will only map to assemblies with the extension `_plass`.

Assembly extensions are the most important of this type of parameter, but the differential expression steps also have a weird paramter they pass in: `contrast`. Target files are generated in `ep_utils/generate_all_targets.py` - have a look if you'd like to see how we build targets to pass into the Snakefile, and look for handling`assembly_extensions` and `contrasts`.


## Workflows

Rules are aggregated into Subworkflows and Workflows that include all the rules necessary to run a larger analysis. For example, the `default` eel pond protocol conducts de novo transcriptome assembly, annotation, and quick differential expression analysis.

For now, workflows are specified in the `ep_utils/pipeline_defaults.yaml`. Each workflow (or subworkflow, which are smaller subsets of the workflows) needs two pieces of information: 

  - `include`: rules to include
  - `targets`: the endpoints of workflows that are not used to generate additional files. These are the "targets" we build and pass to the Snakefile.

