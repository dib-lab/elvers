# Sourmash

sourmash is a command-line tool and Python library for computing MinHash sketches from DNA sequences, comparing them to each other, and plotting the results. This allows you to estimate sequence similarity between even very large data sets quickly and accurately. Please see the mash software and the mash paper [(Ondov et al., 2016)](http://biorxiv.org/content/early/2015/10/26/029827) for background information on how and why MinHash sketches work.

Sourmash is [dib-lab](http://ivory.idyll.org/lab/) software! Please see the [sourmash documentation](https://sourmash.readthedocs.io/en/latest/index.html) for more on sourmash. Sourmash 2.0 is coming soon. In the meantime, please cite [Brown and Irber, 2016](https://joss.theoj.org/papers/3d793c6e7db683bee7c03377a4a7f3c9)

## Quickstart: running sourmash via eelpond:

```
./run_eelpond nema-test sourmash_compute
```
Sourmash computes signatures for both assemblies and kmer-trimmed reads. If you haven't conducted quality and kmer trimming, this workflow will conduct those for you. However, you do need to either 1) have already run an assembly, such that a `fasta` file is sitting in the `eelpond/assembly` directory, 2) Run an assembly at the same time, or 3) pass an assembly in via `assemblyinput`

If you have not already run `./run_eelpond nema-test assemble`:

   2) Run trinity assembly at the same time:
   ```
   ./run_eelpond nema-test assemble sourmash
   ```
   3) OR, Pass an assembly in via `assemblyinput`
   ```
   ./run_eelpond assemblyinput sourmash
   ```
   with an assembly in your `yaml` configfile, e.g.:
   ```
   assemblyinput:
     assembly: rna_testdata/nema.fasta
     gene_trans_map:  rna_testdata/nema.fasta.gene_trans_map
     assembly_extension: '_input'
     ```
    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option. If you have these in your configfile, `eelpond` will automatically assume you want to run the `assemblyinput` rules, but it's nice to specify them in the command anyway :).

## Sourmash Command

On the command line, the command eelpond runs for each file is approximately:
```
sourmash compute --scaled 1000 \
  -k 31 input_file -o output.sig
```

## Customizing Sourmash parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a Sourmash configfile you can modify, run:
```
./run_eelpond sourmash.yaml sourmash --build_config
```
The output should be a small `yaml` configfile that contains:
```
  ####################  sourmash  ####################
sourmash:
  k_size: 31
  scaled: 1000
  extra: ''
```
In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra sourmash parameters,  e.g.:
```
  extra: ' --track-abundance '
```
Override default params by modifying any of these lines, and placing them in the config file you're using to run `eelpond`. Here, we just generated params for `sourmash`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file, e.g. `./run_eelpond my-workflow.yaml full --build_config` and edit parameters there.

## Sourmash eelpond rule

We use a slightly modified version of the [sourmash snakemake wrapper](https://github.com/dib-lab/eelpond/blob/master/rules/sourmash/sourmash-wrapper.py) to run Sourmash via snakemake. 

For snakemake afficionados, see our sourmash rules on [github](https://github.com/dib-lab/eelpond/blob/master/rules/sourmash/sourmash.rule).
