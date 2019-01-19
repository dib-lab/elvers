# Trimmomatic

We use [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.36) to trim off residual Illumina adapters that were left behind after demultiplexing.

Here we use a set TruSeq Illumina adapters. However, if running this on your own data and you know you have different adapters, you'll want to input them in the configfile (see `params` section, below). If you're using the right adapters,  you should see that some of the reads are trimmed; if they’re not, you won’t see anything get trimmed.

See excellent paper on trimming parameters by [MacManes 2014](https://www.frontiersin.org/articles/10.3389/fgene.2014.00013/full).

Based on these recommendations by MacManes 2014, we use this command in this pipeline:

```
TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:25
```
## Customizing the trimming pipeline

The default trimming paramters can be overridden by providing the following in the `.yaml` configuration file:

Be sure to modify the trimming commands as desired

From the [Trimmomatic documentation](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf):

```
Trimmomatic is a fast, multithreaded command line tool that can be used to trim and crop
Illumina (FASTQ) data as well as to remove adapters. These adapters can pose a real problem
depending on the library preparation and downstream application.

There are two major modes of the program: Paired end mode and Single end mode. The
paired end mode will maintain correspondence of read pairs and also use the additional
information contained in paired reads to better find adapter or PCR primer fragments
introduced by the library preparation process.

Trimmomatic works with FASTQ files (using phred + 33 or phred + 64 quality scores,
depending on the Illumina pipeline used). Files compressed using either „gzip‟ or „bzip2‟ are
supported, and are identified by use of „.gz‟ or „.bz2‟ file extensions.
```


## Trimmomatic Params

We run `trimmomatic` via snakemake, but on the commandline, the command would look like this:

```
TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:15 \
    MINLEN:25
```

The trimming command can be modified via the configfile. This is what that same command looks
like in our config:
```
trimmomatic:
  trim_cmd:  "ILLUMINACLIP:{}:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:25"
```
To get the trimmomatic config you can modify, run:
```
./run_eelpond trimmomatic.config --build_config trimmomatic
```
The output should look like this:
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


In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra trimmomatic parameters, e.g.:
```
  extra: '--someflag someparam --someotherflag thatotherparam'
```

## Trimmomatic Snakemake Wrappers

We use a local copies of the [trimmomatic snakemake wrappers](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/trimmomatic.html) to run Trimmomatic.

## eelpond rule

Check out the `eelpond` trimmomatic rule on [github](https://github.com/dib-lab/eelpond/blob/master/rules/trimmomatic/trimmomatic.rule)

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
