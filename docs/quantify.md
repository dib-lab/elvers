# Quantify Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting withing `eelpond`. The "quantify" subworkflow conducts read quality trimming and salmon quantification. It requires an assembly to be provided, either by running an assembly or providing one in your configfile. At the moment, this workflow consists of:
 
  - [get_data](get_data.md) - an `eelpond` utility
  - [trimmomatic](trimmomatic.md)
  - [salmon](salmon.md)


## Quickstart

If you've generated an assembly, even if you've already run `./run_eelpond nema-test assemble`:

   1) "Run" trinity assembly at the same time. If you've already run the assembly, `eelpond` will just locateyour assembly file for `quantify`.

   ```
   ./run_eelpond nema-test assemble quantify
   ```

   2) OR, Pass an assembly in via `assemblyinput` with an assembly in your `yaml` configfile, e.g.: 
   
   ```
   ./run_eelpond assemblyinput quantify
   ```
   
   In the configfile:

    assemblyinput:
      assembly: rna_testdata/nema.fasta
      gene_trans_map:  rna_testdata/nema.fasta.gene_trans_map #optional
      assembly_extension: '_input'

    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have a gene to transcript map, please specify it as well.   If not, delete this line from your `config`. The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Feel free to use `_input`, as   specified above, or pick something equally simple yet more informative. **Note:
    Please don't use additional underscores (`_`) in this extension!**. For more details, see the [assemblyinput documentation](assemblyinput.md). 


## Configuring the quantify subworkflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](about_and_configure.md).

If you want to add the `quantify` program parameters to a previously built configfile, run:
```
./run_eelpond config.yaml quantify --print_params
```

A small set of parameters should print to your console:

```
 ####################  quantify  ####################
get_data:
  download_data: false
  use_ftp: false
trimmomatic:
  adapter_file:
    pe_path: ep_utils/TruSeq3-PE-2.fa
    se_path: ep_utils/TruSeq3-SE.fa
  extra: ''
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25
salmon:
  input_trimmomatic_trimmed: True
    index_params:
      extra: '' 
    quant_params:
      libtype: A
      extra: ''
  #######################################################
```

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](about_and_configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [salmon](salmon.md)