import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
from common.utils import * #get_params,is_se,is_single_end,generate_data_targs,generate_base_targs # *

min_version("5.1.2") #minimum snakemake version

# read in sample info --> use pyflakes instead to check yaml or tsv format?
samples = pd.read_table(config["samples"],dtype=str).set_index(["sample", "unit"], drop=False)
# if units: # add check for unit values (ignore all units if no unit values) 
SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()
validate(samples, schema="schemas/samples_v2.schema.yaml") # new version

# set up dirs, basename
BASE = config.get('basename','eelpond')
experiment_suffix = config.get('experiment_suffix')

if experiment_suffix:
    OUT_DIR = BASE + "_out_" + experiment_suffix
else:
    OUT_DIR = BASE + '_out'

# set up dirs
RULES_DIR = 'rules'
LOGS_DIR = join(OUT_DIR, 'logs')
DATA_DIR = config.get('data_directory', join(OUT_DIR, 'data'))
ANIMALS_DIR = join("common", "animals")
READS_DIR = join(OUT_DIR, 'untrimmed')

TRIM_DIR = join(OUT_DIR, 'trimmed')
QC_DIR = join(OUT_DIR, 'qc')
ASSEMBLY_DIR = join(OUT_DIR, 'assembly')
QUANT_DIR = join(OUT_DIR, 'quant')
KHMER_TRIM_DIR = join(OUT_DIR, 'khmer')
SOURMASH_DIR = join(OUT_DIR,'sourmash')
BUSCO_DIR = join(OUT_DIR,'busco')
DSEQ2_DIR = join(OUT_DIR,'deseq2')
EDGER_DIR = join(OUT_DIR, 'edgeR')
ANNOT_DIR = join(OUT_DIR,'annotation')
BT2_DIR = join(OUT_DIR,'bowtie2')

#print_animal
animal_targs = [join(ANIMALS_DIR,"octopus"),join(ANIMALS_DIR,"fish")]

# check for replicates ** need to change with new samples scheme
# change this replicate check to work with single samples file
replicates = True
num_reps = samples['condition'].value_counts().tolist()
if any(x < 2 for x in num_reps):
    replicates = False

#### snakemake ####

onstart: 
    shell('cat {animal_targs[0]}')
    print('-----------------------------------------------------------------')
    print('Welcome to the Eel Pond, de novo transcriptome assembly pipeline.')
    print('-----------------------------------------------------------------')

onsuccess:
    #if "verbose" in config and config["verbose"]:
    print("\n--- Eel Pond Workflow executed successfully! ---\n")
    shell('cat {animal_targs[1]}')


# download or softlink data
if config.get('download_data', False):
    include: join(RULES_DIR, 'utils', 'ftp.rule')
else:
    include: join(RULES_DIR, 'utils', 'link_data.rule')

    data_ext = [".fq.gz", ".fq.gz"]
    data_targs = generate_data_targs(DATA_DIR, SAMPLES, data_ext)

# program rules
include: join(RULES_DIR,'fastqc', 'fastqc.rule')
include: join(RULES_DIR, 'trimmomatic', 'trimmomatic-pe.rule')
include: join(RULES_DIR, 'utils', 'common.rule')

#diginorm = config.get('diginorm', True)
#include: join(RULES_DIR, 'khmer','khmer_no_diginorm.rule')
include: join(RULES_DIR, 'khmer','khmer.rule')
include: join(RULES_DIR, 'salmon', 'salmon.rule')
include: join(RULES_DIR, 'bowtie2', 'bowtie2.rule')
include: join(RULES_DIR, 'dammit', 'dammit.rule')
include: join(RULES_DIR, 'sourmash', 'sourmash.rule')
include: join(RULES_DIR, 'deseq2', 'deseq2.rule')
include: join(RULES_DIR, 'busco', 'busco.rule')
#include: join(RULES_DIR, 'edgeR', 'edgeR.rule')
#include: join(RULES_DIR, 'sourmash', 'push_sigs.rule')


# embed extensionis within the rules files? --> would be way better. Then just build targs from that.
fastqc_ext =  ['_fastqc.zip','_fastqc.html', '_trimmed_fastqc.zip','_trimmed_fastqc.html']
fastqc_targs = generate_data_targs('qc', SAMPLES, fastqc_ext)

trim_ext = [".trim.fq.gz", ".se.trim.fq.gz"]
trim_targs = generate_data_targs('trimmomatic', SAMPLES, trim_ext)

busco_ext = ['']
busco_targs = generate_base_targs(BUSCO_DIR, 'run_busco_' + BASE, busco_ext)

sourmash_ext = ['.sig'] 
sourmash_targs = generate_base_targs('sourmash', BASE, sourmash_ext)

assembly_targs=[]
assemblies=[]
#if config.get('assembly_program', '').lower() == 'trinity': # enable list of assembly programs?
trinity_ext = ['_trinity.fasta', '_trinity.fasta.gene_trans_map']
trinity_targs = generate_base_targs(ASSEMBLY_DIR, BASE, trinity_ext)
assemblies+=['trinity']
assembly_targs+=trinity_targs

dammit_ext = ['.fasta.dammit.gff3', '.fasta.dammit.fasta']
dammit_targs = generate_base_targs(ANNOT_DIR, BASE, dammit_ext)

salmon_ext = ['/quant.sf', '/lib_format_counts.json']
salmon_targs = generate_data_targs(QUANT_DIR, BASE, salmon_ext, ends = [''])

bt2_ext = ['.bam'] 
bt2_targs = generate_data_targs('bowtie2', BASE, bt2_ext, ends = [''])

deseq2_targs =  ['']#get_targets(units,BASE,DSEQ2_DIR, conf = config)


# eelpond target rules
rule quality_trim:
    input: trim_targs #fastqc_targs + trim_targs

rule assemble:
    input: assembly_targs

rule assembly_quality:
    input: busco_targs + sourmash_targs

rule annotation:
    input: dammit_targs

rule quantification:
    input: salmon_targs

rule bt_map:
    input: bt2_targs

rule diffexpression:
    input: deseq2_targs


#if input_assembly:
#    assert config['assembly_input']['assembly'] is not None, "chosen workflow requires transcriptome assembly as input" 
#    assert config['assembly_input']['gene_trans_map'] is not None, "chosen workflow requires transcriptome gene transcript map as input" 

    # read processing options
	#kmer_trim = config.get('kmer_trim', True)
    #diginorm = config.get('diginorm', True)
    
	#if kmer_trim:
        # add se khmer option back in 
	#	khmer_pe_ext = ['_1.khmer.fq.gz', '_2.khmer.fq.gz', '.paired.khmer.fq.gz', '.single.khmer.fq.gz']
    #    khmer_targs = generate_data_targs(KHMER_TRIM_DIR, SAMPLES, khmer_pe_ext, ends = [""])
        #TARGETS += khmer_targs
    
#if input_assembly:
#    include: 'rules/assemblyinput/assemblyinput.rule'
#    assemblyinput_ext = ['.fasta', '.fasta.gene_trans_map']
#    assemblyinput_targs = generate_base_targs(ASSEMBLY_DIR, BASE, assemb_input_ext)
    #TARGETS += assemblyinput_targs

        #include: 'rules/edgeR/edgeR_no_replicates.rule'
        #from rules.edgeR.edgeR_targets import get_targets
        #edgeR_targs = get_targets(units,BASE,EDGER_DIR, conf = config)
        #TARGETS += edgeR_targs


##### report #####

#report: "report/workflow.rst"
#shell('cat {animal_targs[1]}')

