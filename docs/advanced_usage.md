# Advanced Usage


Each independent step is split into a smaller workflow that can be run independently, if desired, e.g. `elvers examples/nema preprocess`. Individual tools can also be run independently, see [Advanced Usage](advanced_usage.md).

See the help, here:
```
elvers -h
```
**available workflows:**

  - preprocess: Read Quality Trimming and Filtering (fastqc, trimmomatic)
  - kmer_trim: Kmer Trimming and/or Digital Normalization (khmer)
  - assemble: Transcriptome Assembly (trinity)
  - assemblyinput: Specify assembly for downstream steps
  - annotate : Annotate the transcriptome (dammit, sourmash)
  - quantify: Quantify transcripts (salmon)
  - full: preprocess, kmer_trim, assemble, annotate, quantify

You can see the available workflows (and which programs they run) by using the `--print_workflows` flag:
```
elvers examples/nema.yaml --print_workflows
```

*more info coming soon*

