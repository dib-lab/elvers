# Assembly Input Utility Rule


If you want to use the downstream steps such as read quantification and differential expression on a previous assembly, you can use the `assemblyinput` utility to pass this assembly into `elvers`.


To do this, specify the file path for your input assembly (and optionally, a tab-separated gene to transcript map) in the config file.

# Specify Input Assembly

The following will build a config for just the assemblyinput rule. 
```
./run_elvers new_config.yaml assemblyinput --build_config
```

In here, you'll see a section for "assemblyinput" parameters that looks like this:

```
  ####################  assemblyinput  ####################
assemblyinput:
  assembly: examples/nema.assembly.fasta
  gene_trans_map: examples/nema.assembly.fasta.gene_trans_map
  assembly_extension: _input
```

These filename correspond to a test Trinity assembly we provide for the `nema` test data. To use your own files, specify the file path to your `fasta` file. If you have a gene to transcript map, please specify it as well. If not, delete this line from your `config`. The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Feel free to use `_input`, as specified above, or pick something equally simple yet
more informative. **Note: Please don't use additional underscores (`_`) in this extension!**


*Coming soon: option to download this assembly instead*

## Output Files

The output of the `assemblyinput` step is your assembly (and optional gene-transcript map) copied into a subdirectory (`assembly`) within your output directory (BASENAME_out). 


