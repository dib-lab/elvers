import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
from common.utils import get_params

min_version("5.1.2") #minimum snakemake version

# read in & validate sample info 
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")
units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

# build file extensions from suffix info (+ set defaults)
base = config.get('basename','eelpond') 
experiment_suffix =config.get('experiment', '')

# build directory info --> later set all these from config file(s)
#folders = config['directories']

OUT_DIR = '{}_out{}'.format(base, experiment_suffix)
LOGS_DIR = join(OUT_DIR, 'logs')
TRIM_DIR = join(OUT_DIR, 'trimmed')
KHMER_TRIM_DIR = join(OUT_DIR, 'khmer')
QC_DIR = join(OUT_DIR, 'qc')
ASSEMBLY_DIR = join(OUT_DIR, 'assembly')
QUANT_DIR = join(OUT_DIR, 'quant')
SOURMASH_DIR = join(OUT_DIR,'sourmash')

# workflow rules

#fastqc
include: 'rules/common.rule'
include: 'rules/fastqc/fastqc.rule'
from rules.fastqc.fastqc_targets import get_targets
fastqc_targs = get_targets(units, base, QC_DIR)
#trimmomatic
include: 'rules/trimmomatic/trimmomatic.rule'
from rules.trimmomatic.trimmomatic_targets import get_targets
trim_targs = get_targets(units, base, TRIM_DIR)
#trinity
include: 'rules/trinity/trinity.rule'
from rules.trinity.trinity_targets import get_targets
trinity_targs = get_targets(units, base, ASSEMBLY_DIR)
#salmon
include: 'rules/salmon/salmon.rule'
from rules.salmon.salmon_targets import get_targets
salmon_targs = get_targets(units, base, QUANT_DIR)
#khmer
include: 'rules/khmer/khmer.rule'
from rules.khmer.khmer_targets import get_targets
khmer_targs = get_targets(units, base, KHMER_TRIM_DIR)
#sourmash
include: 'rules/sourmash/sourmash.rule'
from rules.sourmash.sourmash_targets import get_targets
sourmash_targs = get_targets(base,SOURMASH_DIR)

TARGETS = fastqc_targs + trim_targs + trinity_targs + salmon_targs + sourmash_targs #+ khmer_targs
print(TARGETS)

rule all:
    input: TARGETS 


##### singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### report #####

report: "report/workflow.rst"
