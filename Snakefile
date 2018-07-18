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
experiment_suffix = config.get('experiment_suffix', '')

# build directory info --> later set all these from config file(s)? or just a defaults file?
#folders = config['directories']

ANIMALS_DIR = "common/animals/"
OUT_DIR = '{}_out{}'.format(base, experiment_suffix)
LOGS_DIR = join(OUT_DIR, 'logs')
READS_DIR = join(OUT_DIR, 'untrimmed')
TRIM_DIR = join(OUT_DIR, 'trimmed')
KHMER_TRIM_DIR = join(OUT_DIR, 'khmer')
QC_DIR = join(OUT_DIR, 'qc')
ASSEMBLY_DIR = join(OUT_DIR, 'assembly')
QUANT_DIR = join(OUT_DIR, 'quant')
SOURMASH_DIR = join(OUT_DIR,'sourmash')
BUSCO_DIR = join(OUT_DIR,'busco')
DSEQ2_DIR = join(OUT_DIR,'deseq2')
EDGER_DIR = join(OUT_DIR, 'edgeR')
ANNOT_DIR = join(OUT_DIR,'annotation')

flow = config.get('workflow', 'full')
read_processing,assembly,assembly_quality,annotation,quantification,diffexp,input_assembly = [False]*7 

if flow == 'full': 
    read_processing = True
    assembly = True
    assembly_quality = True
    annotation = True
    quantification = True
    diffexp = True
elif flow =='assembly':
    read_processing = True
    assembly = True
    quality = True
else: 
    input_assembly = True
    assert config['assembly_input']['assembly'] is not None, "chosen workflow requires transcriptome assembly as input" 
    assert config['assembly_input']['gene_trans_map'] is not None, "chosen workflow requires transcriptome gene transcript map as input" 
    if flow == 'annotation':
        annotation = True
    if flow == 'expression':
        read_processing = True
        quantification = True
        diffexp = True

#print_animal
animal_targs = [ANIMALS_DIR+"octopus",ANIMALS_DIR+"fish"]

# workflow rules
TARGETS = []

include: 'rules/common.rule'

if read_processing:
    #fastqc
    include: 'rules/fastqc/fastqc.rule'
    from rules.fastqc.fastqc_targets import get_targets
    fastqc_targs = get_targets(units, base, QC_DIR)
    TARGETS += fastqc_targs
    #trimmomatic
    include: 'rules/trimmomatic/trimmomatic.rule'
    from rules.trimmomatic.trimmomatic_targets import get_targets
    trim_targs = get_targets(units, base, TRIM_DIR)
    TARGETS += trim_targs

if assembly:
    #khmer
    include: 'rules/khmer/khmer.rule'
    from rules.khmer.khmer_targets import get_targets
    khmer_targs = get_targets(units, base, KHMER_TRIM_DIR)
    TARGETS += khmer_targs
    #trinity
    include: 'rules/trinity/trinity.rule'
    from rules.trinity.trinity_targets import get_targets
    trinity_targs = get_targets(units, base, ASSEMBLY_DIR)
    TARGETS += trinity_targs

if input_assembly:
    include: 'rules/assemblyinput/assemblyinput.rule'
    from rules.assemblyinput.assemblyinput_targets import get_targets
    assemblyinput_targs = get_targets(units, base, ASSEMBLY_DIR)
    TARGETS += assemblyinput_targs

if assembly_quality:
    #busco
    include: 'rules/busco/busco.rule'
    from rules.busco.busco_targets import get_targets
    busco_targs = get_targets(units, base, BUSCO_DIR)
    #TARGETS += busco_targs
    #sourmash
    include: 'rules/sourmash/sourmash.rule'
    from rules.sourmash.sourmash_targets import get_targets
    sourmash_targs = get_targets(base,SOURMASH_DIR)
    TARGETS += sourmash_targs

if annotation:
   #dammit
   include: 'rules/dammit/dammit.rule'
   from rules.dammit.dammit_targets import get_targets
   dammit_targs = get_targets(units, base, ANNOT_DIR)
   TARGETS += dammit_targs

if quantification:
    #salmon
    include: 'rules/salmon/salmon.rule'
    from rules.salmon.salmon_targets import get_targets
    salmon_targs = get_targets(units, base, QUANT_DIR)
    TARGETS += salmon_targs

if diffexp:
    #deseq2
    include: 'rules/deseq2/deseq2.rule'
    from rules.deseq2.deseq2_targets import get_targets
    deseq2_targs = get_targets(units,base,DSEQ2_DIR, conf = config)
    TARGETS += deseq2_targs
    include: 'rules/edgeR/edgeR.rule'
    from rules.edgeR.edgeR_targets import get_targets
    edgeR_targs = get_targets(units,base,EDGER_DIR, conf = config)
    #TARGETS += edgeR_targs

#push_sigs
#include: 'rules/push_sigs.rule'

onstart: 
    shell('cat {animal_targs[0]}')
    print('-----------------------------------------------------------------')
    print('Welcome to the Eel Pond, de novo transcriptome assembly pipeline.')
    print('-----------------------------------------------------------------')

    # for testing: snakemake usually prints targets of each rule for us
    # print('Output files to be generated by the pipeline:')
    # print(TARGETS)

rule all:
    input: TARGETS

##### singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

#shell('cat {animal_targs[1]}')

##### report #####

report: "report/workflow.rst"
shell('cat {animal_targs[1]}')

onsuccess:
    #if "verbose" in config and config["verbose"]:
    print("\n--- Eel Pond Workflow executed successfully! ---\n")
    shell('cat {animal_targs[1]}')

