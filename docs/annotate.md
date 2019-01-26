# Annotate Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting withing `eelpond`. The "annotate" subworkflow annotates a transcriptome. It requires an assembly to be provided, either by running an assembly or providing one in your configfile. 

At the moment, this workflow consists of:
 
  - [dammit](dammit.md)

## Quickstart

If you've generated an assembly, even if you've already run `./run_eelpond examples/nema.yaml assemble`:

   1) "Run" trinity assembly at the same time.If you've already run the assembly, `eelpond` will just locateyour assembly file for `annotate`. 
   
   ```
   ./run_eelpond examples/nema.yaml assemble annotate
   ```
  
   2) OR, Pass an assembly in via `assemblyinput` **with an assembly in your `yaml` configfile, e.g.:** 
   
   ```
   ./run_eelpond assemblyinput annotate
   ```
   
   In the configfile:
    
    assemblyinput:
      assembly: examples/nema.assembly.fasta
      gene_trans_map: examples/nema.assembly.fasta.gene_trans_map #optional
      assembly_extension: '_input'
   
   This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have a gene to transcript map, please specify it as well.   If not, delete this line from your `config`. The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Feel free to use `_input`, as   specified above, or pick something equally simple yet more informative. **Note:
    Please don't use additional underscores (`_`) in this extension!**. For more details, see the [assemblyinput documentation](assemblyinput.md). 


## Output files:

Your main output directory will be determined by your config file: by default it is `BASENAME_out` (you specify BASENAME).

Dammit will output files in the `annotation` subdirectory of this output directory. The annotated fasta file will be `ASSEMBLY.dammit.fasta` and the annotation `gff3` file will be `ASSEMBLY.dammit.gff3`. Dammit will also produce a number of intermediate files that will be contained within an `ASSEMBLY.fasta.dammit` folder.

## Configuring the annotate subworkflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](about_and_configure.md).

If you want to add the `annotate` program parameters to a previously built configfile, run:
```
./run_eelpond config.yaml annotate --print_params
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

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](about_and_configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [dammit](dammit.md)
