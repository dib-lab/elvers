# Workflows and Tools








## Available Workflows

To see the available workflows, run:
```
./run_elvers examples/nema.yaml --print_workflows
```

You should see something like this:
```
####################  Available Eelpond Workflows  ####################

  default:
	get_data
	trimmomatic
	khmer
	trinity
	fastqc
	dammit
	salmon
	sourmash

  protein_assembly:
	get_data
	trimmomatic
	khmer
	plass
	pear
	fastqc
	paladin
	sourmash

  preprocess:
	get_data
	fastqc
	trimmomatic

  kmer_trim:
	get_data
	trimmomatic
	khmer

  assemble:
	get_data
	trimmomatic
	khmer
	trinity

  annotate:
	dammit

  quantify:
	get_data
	trimmomatic
	salmon

  diffexp:
	get_data
	trimmomatic
	salmon
	deseq2

  sourmash_compute:
	get_data
	trimmomatic
	khmer
	sourmash

  plass_assemble:
	get_data
	trimmomatic
	khmer
	plass

  paladin_map:
	get_data
	trimmomatic
	khmer
	plass
	pear
	paladin

  correct_reads:
	get_data
	trimmomatic
	rcorrector

```
## Available Tools

To see the available integrated programs, run:
```
./run_elvers examples/nema.yaml --print_rules
```
You should see something like this:
```
  ####################  Advanced usage: all available ruless  ####################

	get_data
	trimmomatic
	fastqc
	rcorrector
	khmer
	trinity
	salmon
	deseq2
	dammit
	sourmash
    assemblyinput
	plass
	pear
	paladin

```
