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

# to do: disable "unit" if desired
#def is_se(units,sample, unit = '', end = ''):
#    if unit:
#        return pd.isnull(samples.loc[(sample, unit), "fq2"])
#    else:
#        return pd.isnull(samples.loc[(sample), "fq2"]) #any nulls? what do we want to return?

# sample checks
def is_single_end(sample, unit, end = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def generate_data_targs(outdir, samples, extensions = {}):
    targs = []
    pe_ext = extensions.get('pe', None)
    se_ext = extensions.get('se', None)
    se_names = samples.loc[samples["fq2"].isnull(), 'name'].tolist()
    pe_names = samples.loc[samples["fq2"].notnull(), 'name'].tolist()
    print(se_names)
    if se_ext and len(se_names)>0:
        print('se')
        targs+=[join(outdir, name + e) for e in se_ext for name in se_names]
    if pe_ext and len(pe_names) > 0:
        targs+=[join(outdir, name + e) for e in pe_ext for name in pe_names]
    return targs
    
    #SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()
    #for s in SAMPLES:
        #if is_single_end(s, u):
     #   target_list = target_list + [join(outdir, s + e) for e in exts]
    #return target_list

def generate_base_targs(outdir, basename, extensions):
    target_list = []
    target_list = [join(outdir, basename + e) for e in extensions]
    return target_list

def generate_program_targs(configD, samples, basename):
    # given configD from each program, build targets
    outdir = configD['eelpond_dirname']
    exts = configD['extensions']
    targets = []
    if exts.get('read', None): 
        targets+=generate_data_targs(outdir, samples, exts.get('read'))
    if exts.get('base', None):
        targets+=generate_base_targs(outdir, basename, exts.get('base'))
    return targets

def generate_mult_targs(configD, workflow, samples):
    # pass full config, program names. Call generate_targs to build all
    workflows = configD['eelpond_pipeline']
    targs = []
    base = configD['basename']
    if workflows.get(workflow, None):
        target_rules = configD['eelpond_pipeline'][workflow]['targets']
        for r in target_rules:
            targs += generate_program_targs(configD[r], samples, base)
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


