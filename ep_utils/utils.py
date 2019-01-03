import os
import yaml
import collections
import pandas as pd
from os.path import join

# general utilities
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

def update_nested_dict(d, other):
# Can't just update at top level, need to update nested params
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

def generate_targs(outdir, basename, samples, assembly_exts=[''], base_exts = None, read_exts = None):
    base_targets, read_targets = [],[]
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
    for ext in assembly_exts:
        assemblyname = basename + ext
        if base_exts:    
            base_targets += [join(outdir, assemblyname + e) for e in base_exts]
        # handle read targets that contain assembly info
        read_targets = [t.replace("__assembly__", assemblyname) for t in read_targets]
    return base_targets + read_targets

def generate_program_targs(configD, samples, basename, assembly_exts):
    # given configD from each program, build targets
    outdir = configD['outdir']
    exts = configD['extensions']
    if exts.get('assembly_extensions'): # this program is an assembler or only works with specific assemblies
        assembly_exts = exts.get('assembly_extensions') # override generals with rule-specific assembly extensions
    targets = generate_targs(outdir, basename, samples, assembly_exts, exts.get('base', None),exts.get('read'))
    return targets

def generate_mult_targs(configD, workflow, samples):
    # pass full config, program names. Call generate_program_targs to build each
    workflows = configD['eelpond_pipeline']
    targs = []
    base = configD['basename']
    assembly_exts = configD.get('assembly_extensions', [""])
    if workflows.get(workflow, None):
        target_rules = configD['eelpond_pipeline'][workflow]['targets']
        for r in target_rules:
            targs += generate_program_targs(configD[r]['eelpond_params'], samples, base, assembly_exts)
    return targs
            
# don't need this anymore: just building full config via run_eelpond
def get_params(rule_name, rules_dir='rules'):
    rule_paramsfile = os.path.join(rules_dir, rule_name, rule_name + '_params.yaml') 
    rule_params = {}
    with open(rule_paramsfile, 'r') as stream:
        try:
            paramsD = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
        rule_params= paramsD[rule_name]
    return rule_params


