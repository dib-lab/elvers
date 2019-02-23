# Diffexp Subworkflow

Subworkflows combine tools in the right order to facilitate file targeting withing `eelpond`. The "diffexp" subworkflow conducts read quality trimming, salmon quantification and differential expression analysis. It requires an assembly to be provided, either by running an assembly or providing one in your configfile. If you have a gene-to-transcript map, this will also need to be specified via the `assemblyinput` parameter. If not, you'll need to modify the deseq2 parameters so
that the workflow will not expect this file. 

At the moment, this workflow consists of:
 
  - [get_data](get_data.md) - an `eelpond` utility
  - [trimmomatic](trimmomatic.md)
  - [salmon](salmon.md)
  - [deseq2](deseq2.md)

## Quickstart

If you've generated an assembly, even if you've already run `./run_eelpond examples/nema.yaml assemble`:

   1) "Run" trinity assembly at the same time. If you've already run the assembly, `eelpond` will just locate your assembly file for `diffexp`.
   
   ```
   ./run_eelpond examples/nema.yaml assemble diffexp
   ```

   2) OR, Pass an assembly in via `assemblyinput` with an assembly in your `yaml` configfile, e.g.:
   
   ```
   ./run_eelpond assemblyinput diffexp
   ```
  
   In the configfile:

    ```
    assemblyinput:
      assembly: examples/nema.assemblyfasta
      gene_trans_map:  examples/nema.assembly.fasta.gene_trans_map #optional
      assembly_extension: '_input'
    ```
    
    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have a gene to transcript map, please specify it as well. If not, delete this line from your `config`, and be sure to tell `deseq2` not to use this file by setting `gene_trans_map: False` (see config below).
    
    The `assembly_extension` parameter is important: this is what allows us to build assemblies from several different assemblers on the same dataset. Feel free to use `_input`, as specified above, or pick something equally simple yet more informative. **Note: Please don't use additional underscores (`_`) in this extension!**. For more details, see the [assemblyinput documentation](assemblyinput.md). 


## Configuring the diffexp subworkflow 

To set up your sample info and build a configfile, see [Understanding and Configuring Workflows](configure.md).

If you want to add the `diffexp` program parameters to a previously built configfile, run:
```
./run_eelpond config.yaml diffexp --print_params
```

A small set of parameters should print to your console:

```
 ####################  diffexp  ####################
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
    index_params:
      extra: '' 
    quant_params:
      libtype: A
      extra: ''
deseq2:
  gene_trans_map: True
    # contrasts for the deseq2 results method
  contrasts:
    time0-vs-time6:
      - time0
      - time6
  pca:
    labels:
      # columns of sample sheet to use for PCA
      - condition
  #######################################################
```

Override default params for any program by placing these lines in your `yaml` config file, and modifying values as desired. For more details, see [Understanding and Configuring Workflows](configure.md).For more on what parameters are available, see the docs for each specific program or utility rule:

  - [get_data](get_data.md)
  - [trimmomatic](trimmomatic.md)
  - [salmon](salmon.md)
  - [deseq2](deseq2.md)
