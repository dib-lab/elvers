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
from ep_utils.utils import * 
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
dirs = config['eelpond_directories']

OUT_DIR = dirs['out_dir']
LOGS_DIR = dirs['logs']
assembly_dir = dirs['outdirs']['assemble']
animals_dir = dirs['animals']

#get ascii  animals
animal_targs = glob.glob(join(animals_dir, '*')) # get all ascii animals
animalsD = {os.path.basename(x): x for x in animal_targs}
octopus = animalsD['octopus']
fish = animalsD['fish']

#### snakemake ####
# include rule files
includeRules = config['include_rules']
for r in includeRules:
    include: r

onstart: 
    shell('cat {octopus}')
    print('-----------------------------------------------------------------')
    print('Welcome to the Eel Pond, de novo transcriptome assembly pipeline.')
    print('-----------------------------------------------------------------')

onsuccess:
    print("\n--- Eel Pond Workflow executed successfully! ---\n")
    shell('cat {fish}')

# PREPROCESS
#data_targs = generate_targs(config['link_data'], samples, BASE) 
#fastqc_targs = generate_targs(config['fastqc'], samples, BASE)
#trim_targs = generate_targs(config['trimmomatic'], samples, BASE)

# preprocess targeting rule
rule preprocess:
    input: generate_mult_targs(config, 'preprocess', samples)  #generate_targs(config['trimmomatic'], samples, BASE) + generate_targs(config['fastqc'], samples, BASE) #fastqc_targs + trim_targs

## ASSEMBLY RULES

#if config.get('diginorm', True):
#    include: join(RULES_DIR, 'khmer','khmer.rule')
#else:
#    include: join(RULES_DIR, 'khmer','khmer_no_diginorm.rule')
#khmer_targs = generate_targs(config['khmer'], samples, BASE, ends = [""])

#rule kmer_trim:
#    input: generate_targs(config['khmer'], samples, BASE, ends = [""])

#assmeb_input
#if config['assembly_input']['assembly']:
#    include: join(RULES_DIR, 'utils', 'assemblyinput.rule')
#    assemblyinput_targs = generate_targs(config['assembly_input'], samples, BASE)

#assembly_targs=[]
#assemblies=[]
#if config.get('assembly_program', '').lower() == 'trinity': # enable list of assembly programs?
#trinity_targs = generate_targs(config['trinity'], samples, BASE)

#assemblies+=['trinity']
#assembly_targs+=trinity_targs

#rule assemble:
#    input: generate_targs(config['trinity'], samples, BASE)

## ANNOTATION
#dammit_targs = generate_targs(config['dammit'], samples, BASE)
#rule annotation:
#    input: generate_targs(config['dammit'], samples, BASE)

## QUALITY 
#sourmash_targs = generate_targs(config['sourmash'], samples, BASE)
#busco_targs = generate_targs(config['busco'], samples, 'run_busco_' + BASE)
#push_targs = generate_targs(config['push_sigs'], samples, BASE)
#rule assembly_quality:
#    input: busco_targs + sourmash_targs

## MAPPING AND QUANTIFICATION 
#bowtie2_targs = generate_targs(config['bowtie2'], samples, BASE, ends=[''])
#rule bt_map:
#    input: bowtie2_targs

#salmon_targs = generate_targs(config['salmon'], samples, BASE, ends = [''])
#rule quantification:
#    input: salmon_targs

## DIFFEXP

#deseq2_targs = generate_targs(config['deseq2'], samples, BASE)
#rule diffexpression:
#    input: deseq2_targs

#rule full:
#    input: assembly_targs


#include: join(RULES_DIR, 'edgeR', 'edgeR.rule')
#edgeR_targs = generate_targs(config['edger'], samples, BASE)
        #include: 'rules/edgeR/edgeR_no_replicates.rule'
        #from rules.edgeR.edgeR_targets import get_targets
        #edgeR_targs = get_targets(units,BASE,EDGER_DIR, conf = config)
        #TARGETS += edgeR_targs



##### report #####

#report: "report/workflow.rst"
shell('cat {animal_targs[1]}')

