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
            yamlD = yaml.safe_load(stream) #, Loader=yaml.FullLoader)
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

def handle_reference_input(refInput, config):
    extensions= {}
    program_params = config['get_reference'].get('program_params')
    referencefile = program_params.get('reference', None)
    if not program_params.get('download_ref', False):
        if referencefile:
            assert os.path.exists(referencefile), 'Error: cannot find input reference at {}\n'.format(referencefile)
            sys.stderr.write('\tFound input reference at {}\n'.format(referencefile))
            referencefile = os.path.realpath(referencefile)
        else:
            sys.stderr.write("\n\tError: trying to run `get_reference` workflow, but there's no reference file specified in your configfile. Please fix.\n\n")
        # handle the gene_trans_map
        gtmap = program_params.get('gene_trans_map', '')
        if gtmap:
            assert os.path.exists(gtmap), 'Error: cannot find reference gene_trans_map at {}\n'.format(gtmap)
            sys.stderr.write('\tFound input reference gene-transcript map at {}\n'.format(gtmap))
            extensions = {'base': ['.fasta', '.fasta.gene_trans_map']}
            gtmap = os.path.realpath(gtmap)
        else:
            program_params['gene_trans_map'] = ''
            config['no_gene_trans_map']= True
    # grab user-input reference extension
    input_reference_extension = program_params.get('reference_extension', '_input')
    extensions['reference_extensions'] = [input_reference_extension]
    config['get_reference'] = {'program_params': program_params, 'eelpond_params': {'extensions':extensions}}
    return config, input_reference_extension

def handle_samples_input(config, configfile):
    program_params = config['get_data'].get('program_params')
    samples_file = program_params.get('samples', None)
    if samples_file:
        samples_file = os.path.expanduser(samples_file) #handle ~
        if os.path.exists(samples_file) and not os.path.isdir(samples_file):
            samples_file = os.path.realpath(samples_file)
        else:
            paths = [os.getcwd(), os.path.realpath(os.path.dirname(configfile))]
            for p in paths:
                tryfile = os.path.join(p, samples_file)
                if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                    samples_file = os.path.realpath(tryfile)
                    break
        assert os.path.exists(samples_file), f'Error: cannot find input samples tsv at {samples_file}\n'
        config['get_data']['program_params']['samples'] = samples_file
        sys.stderr.write(f'\tFound input samples tsv file at {samples_file}\n')
    else:
        sys.stderr.write("\n\tError: trying to run `get_data` workflow, but the samples tsv file is not specified in your configfile. Please fix.\n\n")
    return config

def generate_targs(outdir, basename, samples, ref_exts=[''], base_exts = None, read_exts = None, other_exts = None, contrasts = []):
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

    # handle base targets, e.g. refname_ext.fasta or refname_ext.stats
    read_targs = []
    base_targs = []
    # build base targets (using reference extensions)
    if ref_exts:
        for ref_e in ref_exts:
            refname = basename + ref_e
            if base_exts:
                for base_extension in base_exts:
                    base_targets += [join(outdir, refname + e) for e in base_exts]
            if read_targets:
                # handle read targets that contain reference info
                read_targs+= [t.replace("__reference__", refname) for t in read_targets] #should return read_targets if nothing to replace
    else:
        read_targs = read_targets
    #handle contrasts within the base targets
    if contrasts and base_targets:
        for c in contrasts:
            base_targs = [t.replace("__contrast__", c) for t in base_targets]
    else:
        base_targs = base_targets

    # handle outputs with no name into (e.g. multiqc)
    if other_exts:
        # no extension, just single filename that we don't want to generate for multiple reference extensions
        other_targs = [join(outdir, e) for e in other_exts]
    return base_targs + read_targs + other_targs

def generate_program_targs(configD, samples, basename,ref_exts, contrasts):
    # given configD from each program, build targets
    outdir = configD['outdir']
    exts = configD['extensions']
    if exts.get('reference_extensions'): # this program is an assembler or only works with specific assemblies
        ref_exts = exts.get('reference_extensions', ['']) # override generals with rule-specific reference extensions
    targets = generate_targs(outdir, basename, samples, ref_exts, exts.get('base', None),exts.get('read'), exts.get('other'), contrasts)
    return targets

def generate_all_targs(configD, samples=None):
    # pass full config, program names. Call generate_program_targs to build each
    workflows = configD['elvers_workflows']
    targs = []
    base = configD['basename']
    ref_exts = configD.get('reference_extensions', [""])
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
        targs += generate_program_targs(configD[r]['elvers_params'], samples, base, ref_exts, contrasts)
    targs = list(set(targs))
    return targs

def get_params(rule_name, rule_dir='rules'):
    # pass in a rule name & the directory that contains its paramsfile.
    # Return paramsD
    rule_paramsfile = os.path.join(rule_dir, 'params.yml')
    rule_params = {}
    with open(rule_paramsfile, 'r') as stream:
        try:
            paramsD = yaml.safe_load(stream) #, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)
        rule_params= paramsD[rule_name] # in case you have multiple programs in single params.yml
    return rule_params
