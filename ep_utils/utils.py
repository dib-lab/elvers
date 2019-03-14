import os
import sys
import yaml
import collections
import pandas as pd
from os.path import join


# general utilities
def find_Snakefile(workdir):
    snakefile = os.path.join(workdir, 'Snakefile')
    assert os.path.exists(snakefile), 'Error: cannot find Snakefile at {}\n'.format(snakefile)
    return snakefile

def find_yaml(workdir, filename, name):
    # find the workflow config file
    workflowfile = None
    if os.path.exists(filename) and not os.path.isdir(filename):
        workflowfile = filename
    else:
        for suffix in ('', '.yaml', '.yml'):
            tryfile = os.path.join(workdir, filename + suffix)
            if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                if name != 'pipeline_defaults':
                    sys.stderr.write('\tFound {} file at {}\n'.format(name, tryfile))
                workflowfile = tryfile
                break
    assert workflowfile, f'Error, cannot find specified {name} file {filename}\n\n\n   Use option "--build_config" to build a default {name} at {filename}.\n'
    return workflowfile

def read_yaml(filename):
    with open(filename, 'r') as stream:
        try:
            yamlD = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return yamlD

def write_yaml(yamlD, paramsfile):
    with open(paramsfile, 'w') as params:
        yaml.dump(yamlD, stream=params, indent=2,  default_flow_style=False)

def import_configfile(config_file, configD=None):
    newConfig = {}
    configDict = read_yaml(config_file)
    for key, val in configDict.items():
        if isinstance(val, dict): # note that this means the only dicts allowed in user configs are for programs.
            newConfig[key] = {'program_params': val}
        else:
            newConfig[key] = val
    if configD:
        update_nested_dict(configD, newConfig)
        return configD
    else:
        return newConfig

def update_nested_dict(d, other):
# Can't just update at top level, need to update nested params
# Note that this only keeps keys that already exist in other
#https://code.i-harness.com/en/q/3154af
    for k, v in other.items():
        if isinstance(v, collections.Mapping):
            d_v = d.get(k)
            if isinstance(d_v, collections.Mapping):
                update_nested_dict(d_v, v)
            else:
                d[k] = v.copy()
        else:
            d[k] = v

# sample checks
def is_single_end(sample, unit, end = '', assembly = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def handle_assemblyinput(assembInput, config):
    extensions= {}
    program_params = config['assemblyinput'].get('program_params')
    assemblyfile = program_params.get('assembly', None)
    if assemblyfile:
        assert os.path.exists(assemblyfile), 'Error: cannot find input assembly at {}\n'.format(assemblyfile)
        sys.stderr.write('\tFound input assembly at {}\n'.format(assemblyfile))
        assemblyfile = os.path.realpath(assemblyfile)
    else:
        sys.stderr.write("\n\tError: trying to run `assemblyinput` workflow, but there's no assembly file specified in your configfile. Please fix.\n\n")
    gtmap = program_params.get('gene_trans_map', '')
    if gtmap:
        assert os.path.exists(gtmap), 'Error: cannot find assembly gene_trans_map at {}\n'.format(gtmap)
        sys.stderr.write('\tFound input assembly gene-transcript map at {}\n'.format(gtmap))
        extensions = {'base': ['.fasta', '.fasta.gene_trans_map']}
        gtmap = os.path.realpath(gtmap)
    else:
        program_params['gene_trans_map'] = ''
        config['no_gene_trans_map']= True
    # grab user-input assembly extension
    input_assembly_extension = program_params.get('assembly_extension', '_input')
    extensions['assembly_extensions'] = [input_assembly_extension]
    config['assemblyinput'] = {'program_params': program_params, 'eelpond_params': {'extensions':extensions}}
    return config, input_assembly_extension

def generate_targs(outdir, basename, samples, assembly_exts=[''], base_exts = None, read_exts = None, other_exts = None, contrasts = []):
    base_targets, read_targets, other_targs = [],[],[]
    # handle read targets
    if read_exts:
        pe_ext = read_exts.get('pe', None)
        se_ext = read_exts.get('se', None)
        combine_units = read_exts.get('combine_units', False)
        if combine_units:
            se_names = samples.loc[samples["fq2"].isnull(), 'sample'].tolist()
            pe_names = samples.loc[samples["fq2"].notnull(), 'sample'].tolist()
        else:
            se_names = samples.loc[samples["fq2"].isnull(), 'name'].tolist()
            pe_names = samples.loc[samples["fq2"].notnull(), 'name'].tolist()
        if se_ext and len(se_names)>0:
            read_targets+=[join(outdir, name + e) for e in se_ext for name in se_names]
        if pe_ext and len(pe_names) > 0:
            read_targets+=[join(outdir, name + e) for e in pe_ext for name in pe_names]
    # handle base targets 
    read_targs = []
    if assembly_exts:
        for ext in assembly_exts:
            assemblyname = basename + ext
            if base_exts:
                base_targets += [join(outdir, assemblyname + e) for e in base_exts]
            # handle read targets that contain assembly info
            read_targs+= [t.replace("__assembly__", assemblyname) for t in read_targets]
    else:
        read_targs = read_targets
    base_targs = []
    if contrasts:
        for c in contrasts:
            base_targs +=[t.replace("__contrast__", c) for t in base_targets]
    else:
        base_targs = base_targets
    if other_exts:
        # no extension, just single filename that we don't want to generate for multiple assembly extensions 
        other_targs = [join(outdir, e) for e in other_exts]
    return base_targs + read_targs + other_targs

def generate_program_targs(configD, samples, basename, assembly_exts, contrasts):
    # given configD from each program, build targets
    outdir = configD['outdir']
    exts = configD['extensions']
    if exts.get('assembly_extensions'): # this program is an assembler or only works with specific assemblies
        assembly_exts = exts.get('assembly_extensions') # override generals with rule-specific assembly extensions
    targets = generate_targs(outdir, basename, samples, assembly_exts, exts.get('base', None),exts.get('read'), exts.get('other'), contrasts)
    return targets

def generate_mult_targs(configD, workflow, samples):
    # pass full config, program names. Call generate_program_targs to build each
    workflows = configD['eelpond_workflows']
    targs = []
    base = configD['basename']
    assembly_exts = configD.get('assembly_extensions', [""])
    # add assertion to make sure workflow exists in config!
    if workflows.get(workflow, None):
        target_rules = configD['eelpond_workflows'][workflow]['targets']
        for r in target_rules:
            contrasts = configD[r]['program_params'].get('contrasts', [])
            targs += generate_program_targs(configD[r]['eelpond_params'], samples, base, assembly_exts, contrasts)
    targs = list(set(targs))
    return targs

# replacement for generate_mult_targs, to enable full workflows!
def generate_all_targs(configD, samples):
    # pass full config, program names. Call generate_program_targs to build each
    workflows = configD['eelpond_workflows']
    targs = []
    base = configD['basename']
    assembly_exts = configD.get('assembly_extensions', [""])
    # add assertion to make sure workflow exists in config!
    #if workflows.get(workflow, None):
    target_rules = []
    for flow,info in workflows.items():
        if info.get('targets', None): # this is a workflow, not a single rule
            target_rules += list(info.get('targets'))
        else:
            target_rules += [flow]
    for r in set(target_rules):
        contrasts = configD[r]['program_params'].get('contrasts', [])
        targs += generate_program_targs(configD[r]['eelpond_params'], samples, base, assembly_exts, contrasts)
    targs = list(set(targs))
    return targs
            
def get_params(rule_name, rule_dir='rules'):
    # pass in a rule name & the directory that contains its paramsfile. 
    # Return paramsD
    rule_paramsfile = os.path.join(rule_dir,rule_name+ '_params.yaml')
    rule_params = {}
    with open(rule_paramsfile, 'r') as stream:
        try:
            paramsD = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
        rule_params= paramsD[rule_name]
    return rule_params
