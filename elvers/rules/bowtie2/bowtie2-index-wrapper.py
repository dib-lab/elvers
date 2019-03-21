__author__ = "N Tessa Pierce"
__copyright__ = "Copyright 2018, N Tessa Pierce"
__email__ = "ntpierce@ucdavis.edu"
__license__ = "MIT"


from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

index = snakemake.output[0]
index_basename = index.split('.1.bt2',1)[0]

shell("bowtie2-build --threads {snakemake.threads} {snakemake.params.extra} {snakemake.input} {index_basename} {log}" ) 
