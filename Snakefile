# Snakemake configuration file for running eelpond pipelines.
#
# see script 'run_eelpond' in this directory for a convenient entry point.
#


import os, sys
from os.path import join
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
from ep_utils.utils import * 
import glob

min_version("5.1.2") # minimum snakemake version

from snakemake.remote import FTP, HTTP
FTP = FTP.RemoteProvider()
HTTP = HTTP.RemoteProvider()

# get sample file
try:
    sample_file = config["samples"]
except KeyError:
    print("cannot find 'samples' entry in config file!", file=sys.stderr)
    sys.exit(-1)

# read in sample file
try:
    samples = pd.read_csv(sample_file, dtype=str, sep='\t')
except:
    print("ERROR reading sample file {}.".format(sample_file), file=sys.stderr)
    raise

try:
    samples = samples.set_index(["sample", "unit"], drop=False)
except:
    print("ERROR interpreting sample file {}.".format(sample_file), file=sys.stderr)
    print("(Is this a tab-delimited file?)", file=sys.stderr)
    raise

#validate(samples, schema="schemas/samples_v2.schema.yaml") # new version
samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)

# note, this function *needs* to be in this file, or added somewhere it can be accessed by all rules
def is_single_end(sample, unit, end = '', assembly = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

# check for replicates ** need to change with new samples scheme
# change this replicate check to work with single samples file
#replicates = True
#num_reps = samples['condition'].value_counts().tolist()
#if any(x < 2 for x in num_reps):
#    replicates = False

# Set up directories 
animals_dir = config['eelpond_directories']['animals'] 
animal_targs = glob.glob(join(animals_dir, '*')) # get all ascii animals
animalsD = {os.path.basename(x): x for x in animal_targs}
octopus = animalsD['octopus']
fish = animalsD['fish']

#### snakemake ####

# include all rule files
includeRules = config['include_rules']
for r in includeRules:
    include: r

onstart: 
    shell('cat {octopus}')
    print('-----------------------------------------------------------------')
    print('Welcome to the Eel Pond, de novo transcriptome assembly pipeline.')
    print('-----------------------------------------------------------------')


documentation_base = "https://dib-lab.github.io/eelpond/"

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")
    
    ## this could be done better, but should work for now :)
    print("  Outputs for all workflow steps:\n") 
    for key, val in config.items():
        if isinstance(val, dict):
            if val.get('eelpond_params', None):
                outdir = val['eelpond_params']['outdir']
                sys.stdout.write(f"\t{key}: {outdir}\n")
                docs = documentation_base + key
                sys.stdout.write(f"\t\t     for explanation of this step, see: {docs} \n\n")
    print("\n")
    ##
    shell('cat {fish}')

rule eelpond:
    input: generate_all_targs(config, samples)


#shell('cat {animal_targs[1]}')

