# Khmer k-mer trimming and (optional) diginorm

Before running transcriptome assembly, we recommend doing some kmer spectral error trimming on your dataset, and if you have lots of reads, also performing digital normalization. This has the effect of reducing the computational cost of assembly without negatively affecting the quality of the assembly. We use [khmer](https://khmer.readthedocs.io/) for both of these tasks.

# Khmer Command

Here's the command as it would look on the command line:
```
(interleave-reads.py sample_1.trim.fq sample_2.trim.fq )| \\
(trim-low-abund.py -V -k 20 -Z 18 -C 2 - -o - -M 4e9 --diginorm --diginorm-coverage=20) | \\
(extract-paired-reads.py --gzip -p sample.paired.gz -s sample.single.gz) > /dev/null
```

Trimmed reads are used as input. `trim-low-abund.py ` trims low-abundance reads to a coverage of 18. Here, we also perform digital normalization to a *k*-mer (*k* = 20) coverage of 20.The output files are the remaining reads, grouped as pairs and singles (orphans). For more on `trim-low-abund`, see this [recipe](https://khmer-recipes.readthedocs.io/en/latest/007-variable-coverage-trimming/),

Finally, since Trinity expects separate `left` and `right` files, we use `split-paired-reads.py` to split the interleaved pairs into two files.
```
split-paired-reads.py sample.paired.gz
```

## Customizing Khmer parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a khmer configfile you can modify, run:
```
./run_eelpond khmer.yaml khmer --build_config
```
The output should be a small `yaml` configfile that contains:
```
  ####################  khmer  ####################
khmer:
  C: 3
  Z: 18
  coverage: 20
  diginorm: true
  extra: ''
  ksize: 20
  memory: 4e9
```
Override default params by modifying any of these lines, and placing them in the config file you're using to run `eelpond`.

## eelpond rule

For snakemake afficionados, here's the basic structure of the khmer eelpond rules. Directories and parameters are specified via the configfile, for more, see the rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/khmer/khmer.rule).

```
rule khmer_pe_diginorm:
    """
    kmer trim and diginorm with khmer
    """
    input: unpack(get_trimmed)
    output: 
        paired=join(khmer_dir,'{sample}_{unit}.paired.khmer.fq.gz'),
        single=join(khmer_dir,'{sample}_{unit}.single.khmer.fq.gz'),
    message:
        """--- khmer trimming of low-abundance kmers and digital normalization ---"""
    params:
        k = khmer_params.get('ksize', 20),
        Z = khmer_params.get('Z', 18), 
        C = khmer_params.get('C', 3), 
        memory = khmer_params.get('memory', 4e9),
        cov = khmer_params.get('coverage', 20),
        extra = khmer_params.get('extra', '')
    threads: 10
    log: join(LOGS_DIR, 'khmer/{sample}_{unit}.pe.diginorm.log')
    conda:  'khmer-env.yaml'
    shell: " (interleave-reads.py {input.r1} {input.r2} ) | "
           " (trim-low-abund.py -V -k {params.k} -Z {params.Z} -C {params.C} - -o - -M {params.memory} "
           " --diginorm --diginorm-coverage={params.cov}) | (extract-paired-reads.py --gzip "
           " -p {output.paired} -s {output.single}) > {log}"

rule khmer_split_paired:
    input: join(khmer_dir,'{sample}_{unit}.paired.khmer.fq.gz'),
    output:
        r1=join(khmer_dir, '{sample}_{unit}_1.khmer.fq.gz'),
        r2=join(khmer_dir, '{sample}_{unit}_2.khmer.fq.gz'),
    threads: 2
    log: join(LOGS_DIR, 'khmer/{sample}_{unit}_split_paired' + BASE + '.log')
    conda: "khmer-env.yaml"
    shell: """
        split-paired-reads.py {input} --gzip -1 {output.r1} -2 {output.r2} >> {log}
        """

```
