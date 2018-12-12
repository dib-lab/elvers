# Snakemake configuration file for running eelpond pipelines.
#
# see script 'run_eelpond' in this directory for a convenient entry point.
#
# Quickstart: `conf/run dory-test full`
#

import os
from os.path import join
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
from common.utils import * 
import glob

min_version("5.1.2") #minimum snakemake version

# read in sample info 
samples = pd.read_table(config["samples"],dtype=str).set_index(["sample", "unit"], drop=False)
# to do:  add check for unit values (ignore all units if no unit values) 
#SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()
validate(samples, schema="schemas/samples_v2.schema.yaml") # new version

# check for replicates ** need to change with new samples scheme
# change this replicate check to work with single samples file
#replicates = True
#num_reps = samples['condition'].value_counts().tolist()
#if any(x < 2 for x in num_reps):
#    replicates = False

# grab dirs, basename built in run_eelpond
BASE = config['basename']
OUT_DIR = config['out_dir']
LOGS_DIR = config['logs_dir']
ASSEMBLY_DIR = config['assembly_dir']

#get ascii  animals
#animal_targs = [join(config['animals_dir'],"octopus"),join(config['animals_dir'],"fish")]
animal_targs = glob.glob(join(config['animals_dir'], '*')) # get all ascii animals
animalsD = {os.path.basename(x): x for x in animal_targs}
octopus = animalsD['octopus']
fish = animalsD['fish']
#print(animalsD.keys())

#### snakemake ####
# include rule files
#rules = glob.glob(join(config['rules_dir'], '*/*.rule'))
#print(rules)
#for r in rules:
#    include: r


onstart: 
    shell('cat {octopus}')
    print('-----------------------------------------------------------------')
    print('Welcome to the Eel Pond, de novo transcriptome assembly pipeline.')
    print('-----------------------------------------------------------------')

onsuccess:
    print("\n--- Eel Pond Workflow executed successfully! ---\n")
    shell('cat {fish}')

### include program rules ##
RULES_DIR = config['rules_dir']
include: join(RULES_DIR, 'utils', 'common.rule') # might not need this anymore?




# PREPROCESS RULES
if config.get('download_data', False):
    include: join(RULES_DIR, 'utils', 'ftp.rule')
else:
    include: join(RULES_DIR, 'utils', 'link_data.rule')
data_targs = generate_targs(config['raw_data'], samples, BASE) 

include: join(RULES_DIR,'fastqc', 'fastqc.rule')
fastqc_targs = generate_targs(config['fastqc'], samples, BASE)

include: join(RULES_DIR, 'trimmomatic', 'trimmomatic-pe.rule')
trim_targs = generate_targs(config['trimmomatic'], samples, BASE)

# preprocess targeting rule
rule preprocess:
    input: trim_targs #fastqc_targs + trim_targs

## ASSEMBLY RULES

if config.get('diginorm', True):
    include: join(RULES_DIR, 'khmer','khmer.rule')
else:
    include: join(RULES_DIR, 'khmer','khmer_no_diginorm.rule')
khmer_targs = generate_targs(config['khmer'], samples, BASE, ends = [""])
	#if kmer_trim:
	#	khmer_pe_ext = ['_1.khmer.fq.gz', '_2.khmer.fq.gz', '.paired.khmer.fq.gz', '.single.khmer.fq.gz']
    #    khmer_targs = generate_data_targs(KHMER_TRIM_DIR, SAMPLES, khmer_pe_ext, ends = [""])
        #TARGETS += khmer_targs


if config['assembly_input']['assembly']:
    include: join(RULES_DIR, 'utils', 'assemblyinput.rule')
    assemblyinput_targs = generate_targs(config['assembly_input'], samples, BASE)

assembly_targs=[]
assemblies=[]
#if config.get('assembly_program', '').lower() == 'trinity': # enable list of assembly programs?
trinity_targs = generate_targs(config['trinity'], samples, BASE)

assemblies+=['trinity']
assembly_targs+=trinity_targs

rule assemble:
    input: assembly_targs

## ANNOTATION RULES

include: join(RULES_DIR, 'dammit', 'dammit.rule')
dammit_targs = generate_targs(config['dammit'], samples, BASE)

rule annotation:
    input: dammit_targs

## QUALITY RULES

include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
sourmash_targs = generate_targs(config['sourmash'], samples, BASE)

include: join(RULES_DIR, 'busco', 'busco.rule')
busco_targs = generate_targs(config['busco'], samples, 'run_busco_' + BASE)

#include: join(RULES_DIR, 'sourmash', 'push_sigs.rule')
#push_targs = generate_targs(config['push_sigs'], samples, BASE)

rule assembly_quality:
    input: busco_targs + sourmash_targs

## MAPPING AND QUANTIFICATION RULES

include: join(RULES_DIR, 'salmon', 'salmon.rule')
salmon_targs = generate_targs(config['salmon'], samples, BASE, ends = [''])

include: join(RULES_DIR, 'bowtie2', 'bowtie2.rule')
bowtie2_targs = generate_targs(config['bowtie2'], samples, BASE, ends=[''])

rule bt_map:
    input: bowtie2_targs

rule quantification:
    input: salmon_targs

## DIFFEXP RULES

include: join(RULES_DIR, 'deseq2', 'deseq2.rule')
deseq2_targs = generate_targs(config['deseq2'], samples, BASE)

#include: join(RULES_DIR, 'edgeR', 'edgeR.rule')
#edgeR_targs = generate_targs(config['edger'], samples, BASE)
        #include: 'rules/edgeR/edgeR_no_replicates.rule'
        #from rules.edgeR.edgeR_targets import get_targets
        #edgeR_targs = get_targets(units,BASE,EDGER_DIR, conf = config)
        #TARGETS += edgeR_targs

rule diffexpression:
    input: deseq2_targs


rule full:
    input: assembly_targs





##### report #####

#report: "report/workflow.rst"
shell('cat {animal_targs[1]}')

