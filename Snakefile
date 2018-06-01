import os
import numpy as np
import pandas as pd
from utils import container_image_is_external, container_image_name
from snakemake.utils import validate, min_version
min_version("5.1.2") #minimum snakemake version

# add validation. see https://github.com/snakemake-workflows/rna-seq-star-deseq2
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["samples"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="schemas/units.schema.yaml")

units['read_type'] = np.where(units['fq2'].isna(), 'se', 'pe') #PE,SE
#paired = units[units.read_type == 'pe']
#single = units[units.read_type == 'se']

def is_single_end(sample, unit, end):
    return units.loc[(sample,unit)]['read_type'] == 'se' 

# build file extensions from suffix info (+ set defaults)
base = config.get('basename','eelpond')
experiment_suffix = config.get('experiment', '')
readfilt = config['read_filtering']
trim_suffix = readfilt.get('trim_suffix', 'trimmed')
read_qual = readfilt['quality_assessment']

# build directory info
OUT_DIR = '{}_out{}'.format(base, experiment_suffix)
LOGS_DIR = os.path.join(OUT_DIR, 'logs')
TRIM_DIR = os.path.join(OUT_DIR, trim_suffix)
QC_DIR = os.path.join(OUT_DIR, 'qc')

# read filtering workflow rules
include: os.path.join('rules',"fastqc.rule")
add_fastqc_targs = True
include: os.path.join("rules","trimmomatic.rule")

#temporary get-targs function. move to separate file.
def get_targets(units, TRIM_DIR, QC_DIR, fastqc_targs = False): 
    """
    Use the sample info provided in the tsv file 
    to generate required targets for each workflow
    """
    #pre_trim_fastqc_targets = []
    post_trim_fastqc_targets = []
    trim_targets = []
    for s, u in units.iterrows():
        sample, unit, read_type = u['sample'],u['unit'],u['read_type']
        trim_targets.append(os.path.join(TRIM_DIR, '{}_{}_1.trim.fq.gz'.format(sample,unit)))
        if fastqc_targs:
            post_trim_fastqc_targets.append(os.path.join(QC_DIR, '{}_{}_1.trim_fastqc.zip'.format(sample,unit)))
            post_trim_fastqc_targets.append(os.path.join(QC_DIR, '{}_{}_1.trim_fastqc.html'.format(sample,unit)))
        if read_type == 'pe':
            trim_targets.append(os.path.join(TRIM_DIR, '{}_{}_2.trim.fq.gz'.format(sample,unit)))
            if fastqc_targs:
                post_trim_fastqc_targets.append(os.path.join(QC_DIR, '{}_{}_2.trim_fastqc.zip'.format(sample,unit)))
                post_trim_fastqc_targets.append(os.path.join(QC_DIR, '{}_{}_2.trim_fastqc.html'.format(sample,unit)))
    targs = trim_targets + post_trim_fastqc_targets
    #return targs
    return post_trim_fastqc_targets



TARGETS = get_targets(units, TRIM_DIR,QC_DIR, add_fastqc_targs)
#print(TARGETS)
rule all:
    input: TARGETS


