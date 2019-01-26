# Paladin_Map Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting within `eelpond`. The "paladin_map" subworkflow conducts read quality trimming and paladin mapping to a protein reference. It requires an assembly to be provided, either by running an assembly or providing one in your configfile. At the moment, this workflow consists of:
 
  - [get_data](get_data.md) - an `eelpond` utility
  - [trimmomatic](trimmomatic.md)
  - [pear](pear.md)
  - [paladin](paladin.md)


## Quickstart

If you've generated an assembly, even if you've already run `./run_eelpond examples/nema.yaml assemble`:

   1) "Run" trinity assembly at the same time. If you've already run the assembly, `eelpond` will just locateyour assembly file for `paladin_map`.
   
   ```
   ./run_eelpond examples/nema.yaml assemble paladin_map
   ```

   2) OR, Pass an assembly in via `assemblyinput`, with an assembly in your `yaml` configfile, e.g.:
   ```
   ./run_eelpond assemblyinput paladin_map
   ```

   In the configfile:

    assemblyinput:
      assembly: examples/nema.assembly.fasta
      gene_trans_map:  examples/nema.assembly.fasta.gene_trans_map #optional
      assembly_extension: '_plass'
    
    
This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have a gene to transcript map, please specify it as well.   If not, delete this line from your `config`. The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Feel free to use `_input`, as   specified above, or pick something equally simple yet more informative. **Note:
    Please don't use additional underscores (`_`) in this extension!**. For more details, see the [assemblyinput documentation](assemblyinput.md). 

## Configuring the paladin_map subworkflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](about_and_configure.md).

If you want to add the `paladin_map` program parameters to a previously built configfile, run:
```
./run_eelpond config.yaml paladin_map --print_params
```

A small set of parameters should print to your console:

```
 ####################  paladin_map  ####################
get_data:
  download_data: false
  use_ftp: false
trimmomatic:
  adapter_file:
    pe_path: ep_utils/TruSeq3-PE-2.fa
    se_path: ep_utils/TruSeq3-SE.fa
  extra: ''
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25
pear:
  input_kmer_trimmed: false
  input_trimmomatic_trimmed: true
  max_memory: 4G
  pval: 0.01
  extra: ''
paladin:
  alignment_params:
    extra: ''
    f: 125
  index_params:
    reference_type: '3'
    gff_file: ''
  #######################################################
```

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](about_and_configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [pear](pear.md)
  - [paladin](paladin_map.md)
