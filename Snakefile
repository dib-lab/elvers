import os
from os.path import join
import yaml
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version
from common.utils import get_params

min_version("5.1.2") #minimum snakemake version

# read in sample info
samples = pd.read_table(config["samples"],dtype=str).set_index(["sample", "unit"], drop=False)
SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()
validate(samples, schema="schemas/samples_v2.schema.yaml") # new version

#util functions
def is_single_end(sample, unit, end = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def generate_data_targs(outdir, samples, extensions, ends = ["_1", "_2"]):
    target_list = []
    # to do: add paired vs single end check here to generate `ends`
    exts = [x+y for x in ends for y in extensions]
    for s in samples:
        #if is_single_end(s, u):
        target_list = target_list + [join(outdir, s + e) for e in exts]
    return target_list

def generate_base_targs(outdir, basename, extensions):
    target_list = []
    target_list = [join(outdir, basename + e) for e in extensions]
    return target_list

# set up dirs, basename
BASE = config.get('basename','eelpond')
experiment_suffix = config.get('experiment_suffix')

if experiment_suffix:
    OUT_DIR = BASE + "_out_" + experiment_suffix
else:
    OUT_DIR = BASE + '_out'

# check for replicates ** need to change with new samples scheme
replicates = True
num_reps = samples['condition'].value_counts().tolist()
if any(x < 2 for x in num_reps):
    replicates = False

#dirs = config['directories']
RULES_DIR = 'rules'
LOGS_DIR = join(OUT_DIR, 'logs')
DATA_DIR = config.get('data_directory', join(OUT_DIR, 'data'))
ANIMALS_DIR = "common/animals/"
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

# determine workflow
flow = config.get('workflow', 'full')
read_processing,assembly,assembly_quality,annotation,quantification,diffexp,input_assembly,bt2_map = [False]*8 


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
    if flow == 'bowtie2':
        read_processing = True
        bt2_map = True

#print_animal
animal_targs = [ANIMALS_DIR+"octopus",ANIMALS_DIR+"fish"]

# workflow rules
TARGETS = []

include: 'rules/common.rule'

if read_processing:
    #fastqc
    include: join(RULES_DIR,'fastqc/fastqc.rule')
    fastqc_ext =  ['_fastqc.zip','_fastqc.html', '_trimmed_fastqc.zip','_trimmed_fastqc.html']
    fastqc_targs = generate_data_targs(QC_DIR, SAMPLES, fastqc_ext)
    #TARGETS += fastqc_targs
    
	#trimmomatic
	include: join(RULES_DIR, 'trimmomatic', 'trimmomatic.rule')
    trim_ext = [".trim.fq.gz", ".se.trim.fq.gz"]
    trim_targs = generate_data_targs(TRIM_DIR, SAMPLES, trim_ext)
    #TARGETS += trim_targs

if assembly:
    assemblies=[]
    # read processing options
	kmer_trim = config.get('kmer_trim', True)
    diginorm = config.get('diginorm', True)
    
	if kmer_trim:
        #khmer
		if not diginorm:
		    include: join(RULES_DIR, 'khmer','khmer_no_diginorm.rule'
	    else:
            include: join(RULES_DIR, 'khmer','khmer.rule')
        # add se khmer option back in 
		khmer_pe_ext = ['_1.khmer.fq.gz', '_2.khmer.fq.gz', '.paired.khmer.fq.gz', '.single.khmer.fq.gz']
        khmer_targs = generate_data_targs(KHMER_TRIM_DIR, SAMPLES, khmer_pe_ext, ends = [""])
        #TARGETS += khmer_targs
    
    #trinity
	if config.get('assembly_program', '').lower() == 'trinity': # enable list of assembly programs?
    
        include: join(RULES_DIR, 'trinity', 'trinity.rule')
        trinity_ext = ['_trinity.fasta', '_trinity.fasta.gene_trans_map']
        trinity_targs = generate_base_targs(ASSEMBLY_DIR, BASE, trinity_ext)
        assemblies+=['trinity']
        #TARGETS += trinity_targs

if input_assembly:
    include: 'rules/assemblyinput/assemblyinput.rule'
    assemblyinput_ext = ['.fasta', '.fasta.gene_trans_map']
    assemblyinput_targs = generate_base_targs(ASSEMBLY_DIR, BASE, assemb_input_ext)
    #TARGETS += assemblyinput_targs

if assembly_quality:
    #busco
    include: join(RULES_DIR, 'busco', 'busco.rule'
    busco_ext = ['']
    busco_targs = generate_base_targs(BUSCO_DIR, 'run_busco_' + BASE, busco_ext)
    #TARGETS += busco_targs
    #sourmash
    include: 'rules/sourmash/sourmash.rule'
    sourmash_ext = ['.sig'] 
    sourmash_targs = generate_base_targs(SOURMASH_DIR, BASE, sourmash_ext)
    #TARGETS += sourmash_targs

if annotation:
   #dammit
   include: join(RULES_DIR, 'dammit/dammit.rule')
   dammit_ext = ['.fasta.dammit.gff3', '.fasta.dammit.fasta']
   #from rules.dammit.dammit_targets import get_targets
   #dammit_targs = get_targets(units, BASE, ANNOT_DIR)
   TARGETS += dammit_targs

if quantification:
    #salmon
    include: join(RULES_DIR, 'salmon', 'salmon.rule'
    #from rules.salmon.salmon_targets import get_targets
    salmon_ext = ['/quant.sf', '/lib_format_counts.json']
    salmon_targs = generate_data_targs(QUANT_DIR, BASE, salmon_ext, ends = [''])
    #TARGETS += salmon_targs

if bt2_map:
    #bowtie2
    include: join(RULES_DIR, 'bowtie2', 'bowtie2.rule'
    bt2_ext = ['.bam'] 
    #from rules.bowtie2.bowtie2_targets import get_targets
    bt2_targs = generate_data_targs(BT2_DIR, BASE, bt2_ext, ends = [''])
    #TARGETS += bt2_targs

if diffexp:
    if replicates:
        #deseq2
        include: 'rules/deseq2/deseq2.rule'
        from rules.deseq2.deseq2_targets import get_targets
        deseq2_targs = get_targets(units,BASE,DSEQ2_DIR, conf = config)
        TARGETS += deseq2_targs
        #include: 'rules/edgeR/edgeR.rule'
        #from rules.edgeR.edgeR_targets import get_targets
        #edgeR_targs = get_targets(units,BASE,EDGER_DIR, conf = config)
        #TARGETS += edgeR_targs
    #else:
        #include: 'rules/edgeR/edgeR_no_replicates.rule'
        #from rules.edgeR.edgeR_targets import get_targets
        #edgeR_targs = get_targets(units,BASE,EDGER_DIR, conf = config)
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

