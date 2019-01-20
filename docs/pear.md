# Merging Read Pairs with PEAR

In order to map PE nucleotide reads to protein assemblies with `Paladin`, which cannot yet use paired end reads, we use `PEAR` to merge PE reads prior to mapping.

From the `PEAR` [documentation](https://cme.h-its.org/exelixis/web/software/pear/doc.html):

PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger. It is fully parallelized and can run with as low as just a few kilobytes of memory.

PEAR evaluates all possible paired-end read overlaps and without requiring the target fragment size as input. In addition, it implements a statistical test for minimizing false-positive results. Together with a highly optimized implementation, it can merge millions of paired end reads within a couple of minutes on a standard desktop computer.

## Quickstart: Running PEAR with eelpond

We recommend you run `pear` as part of the `paladin_map` pipeline

```
./run_eelpond nema-test paladin_map
```
This will run trimmomatic trimming prior to PEAR merging of paired end reads and paladin mapping. If you'd like to just run `PEAR`, you can run it by including the `get_data` and `trimmomatic` rules in your pipeline:

```
./run_eelpond nema-test get_data trimmomatic pear
```
Note, `PEAR` only works on paired end reads, as there's nothing to do for single end reads!

## PEAR Command

On the command line, the command eelpond runs for each set of files is approximately:
```
pear -f forward_reads -r reverse_reads \ 
  -p p_value -j snakemake.threads -y max_memory \  
  extra -o output_basename 
```

## Customizing PEAR parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a PEAR configfile you can modify, run:
```
./run_eelpond pear.yaml pear --build_config
```
The output should be a small `yaml` configfile that contains:
```
  ####################  pear  ####################
pear:
  input_kmer_trimmed: false
  input_trimmomatic_trimmed: true
  max_memory: 4G
  pval: 0.01
  extra: ''
```
**Please modify the `max_memory` parameter to fit the needs of your reads and system.**

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass in additional parameters,  e.g.:

```
  extra: ' --some_param that_param '
```
Please see the [PEAR documentation](https://cme.h-its.org/exelixis/web/software/pear/doc.html) for info on the params.

Override default params by modifying any of these lines, and placing them in the config file you're using to run `eelpond`. Here, we just generated params for `PEAR`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file. 

For example, run the following to build default parameters for `trimmomatic`, `fastqc`, `khmer`, `plass`, `pear`, and `paladin`, which you can then edit in `my-workflow.yaml`:
```
./run_eelpond my-workflow.yaml plass_assemble paladin_map --build_config
``` 

## PEAR eelpond rule

We wrote a new [PEAR snakemake wrapper](https://github.com/dib-lab/eelpond/blob/master/rules/pear/pear-wrapper.py) to run PEAR via snakemake. This wrapper has not been added to the official repo yet, but feel free to use as needed.

For snakemake afficionados, see our pear rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/pear/pear.rule).
