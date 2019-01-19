# Trimmomatic

## Quickstart: Running Trimmomatic with eelpond

```
./run_eelpond nema-test preprocess
```
This will run trimmomatic trimming and fastqc on pre-trim and post-trim data.


## Trimmomatic

We use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.36) to trim off residual Illumina adapters that were left behind after demultiplexing.

Here we use a set TruSeq Illumina adapters. However, if running this on your own data and you know you have different adapters, you'll want to input them in the configfile (see `params` section, below). If you're using the right adapters,  you should see that some of the reads are trimmed; if they’re not, you won’t see anything get trimmed.

See excellent paper on trimming parameters by [MacManes 2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full).

## Trimmomatic Command

Based on these recommendations by MacManes 2014, we use this command in this pipeline:

```
TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:15 \
    MINLEN:25
```

However, the trimming command can be extensively modified via the configfile. Here's how the parameters for the command above look in our config:

```
trimmomatic:
  trim_cmd:  "ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25"
  extra: ''
```

## Customizing the trimming parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a trimmomatic config you can modify, run:

```
./run_eelpond trimmomatic.yaml trimmomatic --build_config
```

The output should be a small `yaml` configfile that contains:

```
####################  trimmomatic  ####################
trimmomatic:
  adapter_file:
    pe_name: TruSeq3-PE.fa
    pe_url: https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa
    se_name: TruSeq3-SE.fa
    se_url: https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa
  trim_cmd: ILLUMINACLIP:{}:2:40:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
  extra: ''
```


Override default params by modifying these lines. In addition to changing parameters we've specifically enabled, 
you can modify the "extra" param to pass any extra trimmomatic parameters, e.g.:

```
  extra: '--someflag someparam --someotherflag thatotherparam'
```

Or in Trimmomatic params:

```
  extra: 'HEADCROP:5' # to remove the first 5 bases at the front of the read.
```

Be sure the modified lines go into the config file you're using to run `eelpond`. For more on what parameters are available, see the [Trimmomatic documentation](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

## Trimmomatic Rule 

We use a local copies of the [trimmomatic snakemake wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/trimmomatic.html) to run Trimmomatic.

For snakemake afficionados, here's the basic structure of the trimmomatic eelpond rules. Directories and parameters are specified via the configfile (see the rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/trimmomatic/trimmomatic.rule)).

```
def get_pretrim(w):
    readsD = {}
    if not is_single_end(**w):
        readsD['r1'] = expand(join(data_dir, '{sample}_{unit}_1.fq.gz'),**w)
        readsD['r2'] = expand(join(data_dir, '{sample}_{unit}_2.fq.gz'),**w)
        return readsD
    return expand(join(data_dir, '{sample}_{unit}_1.fq.gz'), **w)

rule trimmomatic_pe:
    """
    Trim reads from the sequencer by trimming or dropping low-quality reads.
    """
    input:
        unpack(get_pretrim)
    output:
        r1=join(trim_dir, "{sample}_{unit}_1.trim.fq.gz"),
        r2=join(trim_dir, "{sample}_{unit}_2.trim.fq.gz"),
        r1_unpaired=join(trim_dir, "{sample}_{unit}_1.se.trim.fq.gz"),
        r2_unpaired=join(trim_dir, "{sample}_{unit}_2.se.trim.fq.gz"),
    message:
        """--- Quality trimming PE read data with Trimmomatic."""
    threads: trim_params.get('cpu', 16)
    params:
        trimmer = (trim_params['trim_cmd'].format(trim_params['adapter_file']['pe_name'])).split(' '),
        extra = ''
    log: join(LOGS_DIR, 'trimmomatic/{sample}_{unit}_pe.log')
    wrapper:
        '0.27.1/bio/trimmomatic/pe'

rule trimmomatic_se:
    """
    Trim reads from the sequencer by trimming or dropping low-quality reads.
    """
    input:
        get_pretrim
    output:
        r1=join(trim_dir, "{sample}_{unit}.se.trim.fq.gz"), 
    message:
        """--- Quality trimming SE read data with Trimmomatic."""
    threads: trim_params.get('cpu', 16)
    params:
        trimmer = (trim_params['trim_cmd'].format(trim_params['adapter_file']['se_name'])).split(' '),
        extra = ''
    log:
        join(LOGS_DIR, 'trimmomatic/{sample}_{unit}_se.log')
    wrapper:
        '0.27.1/bio/trimmomatic/se'
```
