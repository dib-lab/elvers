import os
import yaml
import pandas as pd
from os.path import join

#from common.default_vars import *

def is_se(units,sample, unit = '', end = ''):
    if unit:
        return pd.isnull(units.loc[(sample, unit), "fq2"])
    else:
        return pd.isnull(units.loc[(sample), "fq2"]) #any nulls? what do we want to return?

#util functions
def is_single_end(sample, unit, end = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def generate_data_targs(outdir, samples, extensions, ends = ["_1", "_2"]):
    target_list = []
    # to do: add paired vs single end check here to generate `ends`
    exts = [x+y for x in ends for y in extensions]
    SAMPLES = (samples['sample'] + '_' + samples['unit']).tolist()
    for s in SAMPLES:
        #if is_single_end(s, u):
        target_list = target_list + [join(outdir, s + e) for e in exts]
    return target_list

def generate_base_targs(outdir, basename, extensions):
    target_list = []
    target_list = [join(outdir, basename + e) for e in extensions]
    return target_list

def generate_targs(configD, samples, basename, ends = ["_1", "_2"]):
    # given configD from each program, build targets
    outdir = configD['eelpond_dirname']
    exts = configD['extensions']
    targets = []
    if exts.get('read', None):
        targets+=generate_data_targs(outdir, samples, exts.get('read'), ends)
    if exts.get('base', None):
        targets+=generate_base_targs(outdir, basename, exts.get('base'))
    return targets


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


