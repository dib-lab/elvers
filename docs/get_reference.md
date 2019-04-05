# Get Reference Utility Rule

For `elvers` workflows that start with or otherwise utilize a reference, the keyword `get_reference` can be used to provide that reference (and for transcriptomes, an optional tab-separated gene-to-transcript map). If you provide the `get_reference` keyword and information in your configuration file, the `get_reference` utility rule with either download or softlink your reference into your `elvers` directory so it can be utilized for any specified workflows.

# Specify Input reference

The default get_reference parameters are as follows:

```
get_reference:
  reference: examples/nema.assembly.fasta
  gene_trans_map: examples/nema.assembly.fasta.gene_trans_map
  reference_extension: _input
  download_ref: false # download the reference using http (or ftp)
  use_ftp: false # switch download method from http to ftp
```
( to see these, you can run `elvers config get_reference --print_params`)

The default filenames correspond to a test Trinity transcriptome we provide for the `nema` test data. To use your own files, specify the file path to your `fasta` file via the `reference` parameter. If you have a gene to transcript mapi (transcriptomes only), please specify it as well. If not, do not include the `gene_trans_map` line in your `config`. The `reference_extension` parameter is an optional parameter that allows us to use multiple references or assemblies. For example, you can
build a _de novo_ transcriptome with Trinity and input a reference transcriptome and do all downstream steps (e.g. quantification and differential expression) on both assemblies. All assemblies generated via `elvers` rules will have an extension corresponding to the program they were generated with. If you choose to add a `reference_extension`, feel free to use `_input`, as specified above, or pick something equally simple yet
more informative, but **please don't use additional underscores (`_`) in this extension!**


## Output Files

The output of the `get_reference` step is your reference (and optional gene-transcript map) copied into a subdirectory (`reference`) within your output directory (BASENAME_out). 
