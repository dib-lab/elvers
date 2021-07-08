import os
import sys
import yaml
import collections
import pandas as pd
from os.path import join
import glob
from snakemake.utils import validate

# general utilities
def find_Snakefile(workdir):
    snakefile = os.path.join(workdir, 'Snakefile')
    assert os.path.exists(snakefile), 'Error: cannot find Snakefile at {}\n'.format(snakefile)
    return snakefile

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

def read_samples(config):
    elvers_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        samples_file = config["get_data"]["program_params"]["samples"]
    except KeyError:
        print("cannot find 'samples' entry in config file! Please provide a samples file for the `get_data` utility", file=sys.stderr)
        sys.exit(-1)
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            samples = pd.read_csv(samples_file, dtype=str, sep=separator).set_index(["sample", "unit"], drop=False)
            validate(samples, schema=os.path.join(elvers_dir,"schemas/samples_v2.schema.yaml"))
            samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str).set_index(["sample", "unit"], drop=False)
            validate(samples, schema="schemas/samples_v2.schema.yaml") # new version
            samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    # check for single-unit case
    if not (samples['sample'].value_counts() > 1).any():
        config['ignore_units'] = True
    # column 4 is "condition", but can change name
    if (samples.iloc[:, 4].value_counts() < 2).any():
        config['all_replicated'] = False
    return samples, config

# sample checks
def is_single_end(sample, unit, end = '', assembly = ''):
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def find_input_file(filename, name="input file", add_paths=[], add_suffixes = ['.yaml', '.yml'], verbose = False):
    # for any file specified via command line, check if it exists at the current path, if not, try some other paths before returning a helpful error
    found_file = None
    filename = os.path.expanduser(filename) # handle ~!
    paths_to_try = ['', os.getcwd(), os.path.dirname(os.path.abspath(__file__))] + add_paths
    suffixes_to_try = [''] + add_suffixes

    if os.path.exists(filename) and not os.path.isdir(filename):
        found_file = os.path.realpath(filename)
    else:
        for p in paths_to_try:
            for s in suffixes_to_try:
                tryfile = os.path.join(p, filename+ s)
                if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                    found_file = os.path.realpath(tryfile)
                    break
    config_help = ""
    if "config" in name:
        config_help = f"   Use option '--build_config' to build a default {name} at {filename}.\n"

    assert found_file, f'Error, cannot find specified {name} file {filename}\n\n\n' + config_help
    if verbose:
        sys.stderr.write(f'\tFound {name} at {found_file}\n')
    return found_file

def handle_reference_input(config, configfile):
    extensions= {}
    program_params = config['get_reference'].get('program_params')
    referencefile = program_params.get('reference', None)
    if not program_params.get('download_ref', False):
        if referencefile:
            referencefile = find_input_file(referencefile, name="input reference", add_paths = [os.path.realpath(os.path.dirname(configfile))], add_suffixes = ['.fa', '.fasta'])
            program_params['reference'] = referencefile
        else:
            sys.stderr.write("\n\tError: trying to run `get_reference` workflow, but there's no reference file specified in your configfile. Please fix.\n\n")
        # handle the gene_trans_map
        gtmap = program_params.get('gene_trans_map', '')
        if gtmap:
            gtmap = find_input_file(gtmap,"input reference gene_trans_map", add_paths = [os.path.realpath(os.path.dirname(configfile))], add_suffixes = [''])
            program_params['gene_trans_map'] = gtmap
            extensions = {'base': ['.fasta', '.fasta.gene_trans_map']}
        else:
            program_params['gene_trans_map'] = ''
            config['no_gene_trans_map']= True
            extensions = {'base': ['.fasta']}
    # grab user-input reference extension
    input_reference_extension = program_params.get('reference_extension', '')
    extensions['reference_extensions'] = [input_reference_extension]
    config['get_reference'] = {'program_params': program_params, 'elvers_params': {'outputs': {'extensions':extensions}}}
    return config, input_reference_extension

def handle_samples_input(config, configfile):
    program_params = config['get_data'].get('program_params')
    samples_file = program_params.get('samples', None)
    if samples_file:
        samples_file = find_input_file(samples_file, "samples tsv", add_paths = [os.path.realpath(os.path.dirname(configfile))], add_suffixes = ['.tsv'])
        config['get_data']['program_params']['samples'] = samples_file
    else:
        sys.stderr.write("\n\tError: trying to run `get_data` workflow, but the samples tsv file is not specified in your configfile. Please fix.\n\n")
    return config

#def check_workflow(config):
    # This is way too naive. Need to come up with a better version.
    # Maybe we manually check the generated snakemake targs?
    # Or just catch the snakemake error and print better help for which rule needs to be included?.
#    inputs, outputs = [],[]
#    for key, val in config.items():
#        if isinstance(val, dict):
#            if val.get('elvers_params'):
#                inputs += val['elvers_params']['inputs'].get('read', [])
#                inputs += val['elvers_params']['inputs'].get('reference', [])
#                outputs += val['elvers_params']['outputs'].get('read', [])
#                outputs += val['elvers_params']['outputs'].get('reference', [])
                # leaving out "other" inputs/outputs, bc they're never used as inputs (so far)
#    try:
#        set(inputs) == set(outputs)
#    except:
        #well this is uninformative
#        sys.stderr.write("chosen workflow is not valid")


def select_outputs(config):
    # Here we go from "output_options" --> outputs. Mostly, we want to map inputs (program_params) directly to outputs
    # if more than one input, then more than one output. Some of these will have different folders. How do we handle this?
    # might need to get rid of outdir, and just do os.path.dirname on each output? Not sure yet. Or maybe we just select outputs
    # so that we can build all required directories and generate the output targets for snakemake, then put outputs: named_output
    # = same structure as "output_options", just different name, bc it signifies that we're making ALL these outputs.
    reference_extensions = config.get('reference_extensions', [])
    for key, val in config.items():
        if isinstance(val, dict):
            if val.get('elvers_params'):
                if key not in ['get_data', 'get_reference']:
                    # this is a program! proceed
                    inputs = val['program_params']['inputs']
                    outputs = {}
                    ref_exts = []
                    for output_name, output_info in val['elvers_params']['output_options'].items():
                        input_for_this_output = output_info['input']
                        if any([input_for_this_output in inputs, input_for_this_output == 'any']): # choose the output that corresponds to the input going in
                            outputs[output_name] = output_info
                            ref_exts = output_info['extensions'].get('reference_extensions', [])
                            reference_extensions+=ref_exts
                    val['elvers_params']['outputs'] = outputs
                else:
                    val['elvers_params']['outputs'] = val['elvers_params']['output_options']
                del val['elvers_params']['output_options']
                config[key] = val
            # maybe we also want to get rid of the output_options, to minimize clutter
    config['reference_extensions'] = reference_extensions
    return config


def generate_targs(outdir, basename, samples, ref_exts=[''], base_exts = None, read_exts = None, other_exts = None, contrasts = [], ignore_units=False):
    base_targets, read_targets, other_targs = [],[],[]
    # handle read targets
    if read_exts:
        pe_ext = read_exts.get('pe', None)
        se_ext = read_exts.get('se', None)
        combine_units = read_exts.get('combine_units', False)
        if combine_units and not ignore_units:
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
        references = [basename + e for e in ref_exts]
    else:
        references = [basename]
    for refname in references:
        if base_exts:
            base_targets += [join(outdir, refname + e) for e in base_exts]
        if read_targets:
            # handle read targets that contain reference info
            read_targs+= [t.replace("__reference__", refname) for t in read_targets] #should return read_targets if nothing to replace
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

def generate_rule_targs(home_outdir, basename, ref_exts, rule_config, rulename, samples, all_input_options, ignore_units):
    contrasts = rule_config['program_params'].get('contrasts', [])
    # handle outputs (sometimes multiple outputs)
    output_files = []
    outdir = ""
    for outname, output_info in rule_config['elvers_params']['outputs'].items():
        outdir = output_info['outdir']
        out_exts = output_info['extensions']
        if out_exts.get('reference_extensions'): # this program is an assembler or only works with specific assemblies
            out_ref_exts = out_exts.get('reference_extensions', ['']) # override generals with rule-specific reference extensions
        else:
            out_ref_exts = ref_exts

        outputs = generate_targs(outdir, basename, samples, out_ref_exts, out_exts.get('base', None),out_exts.get('read'), out_exts.get('other'), contrasts, ignore_units)
        output_files += outputs
    rule_config['elvers_params']['outputs']['output_files'] = output_files
    rule_config['elvers_params']['outputs']['outdir'] = outdir

    # handle input options!
    if rulename == 'get_data':
        samples_file = rule_config['program_params']['samples'] # should be present (validated) prior to here
        rule_config['elvers_params']['input_options'] = {'get_data': {'indir': os.path.dirname(samples_file), 'input_files': [samples_file]}}
    elif rulename == 'get_reference':
        reference = rule_config['program_params']['reference'] # should be present (validated) prior to here
        rule_config['elvers_params']['input_options'] = {'get_ref': {'indir': os.path.dirname(reference), 'input_files': [reference]}}
    else:
        all_input_exts = {}
        input_options = rule_config['elvers_params'].get('input_options', None)# read, reference, other
        not_found = []
        for input_type, options in input_options.items():
            for option in options:
                try:
                    info = all_input_options[input_type].get(option)
                    all_input_exts[option] = info
                    indir = os.path.join(home_outdir, info.get('outdir', 'input_data'))
                    in_exts = info['extensions']
                    if in_exts.get('reference_extensions'): # this program is an assembler or only works with specific assemblies
                        in_ref_exts = in_exts.get('reference_extensions', ['']) # override generals with rule-specific reference extensions
                    else:
                        in_ref_exts = ref_exts
                    input_files = generate_targs(indir, basename, samples, in_ref_exts, in_exts.get('base', None),in_exts.get('read'), in_exts.get('other'), contrasts)
                    all_input_exts[option]['input_files'] = generate_targs(indir, basename, samples, in_ref_exts, in_exts.get('base', None),in_exts.get('read'), in_exts.get('other'), contrasts)
                except:
                    not_found.append(option)
                if not_found == options:
                    option_list = "\n  " + "\n  ".join(options)
                    sys.stderr.write(f"cannot find input files for {rulename}. Please add a target that produces any of the following: {option_list}")
                    sys.exit(-1)
        rule_config['elvers_params']['input_options'] = all_input_exts
    return rule_config

def generate_inputs_outputs(config, samples=None):
    all_inputs, all_outputs = [],[]
    base = config['basename']
    ref_exts = config.get('reference_extensions', [""])
    home_outdir = config['elvers_directories']['out_dir']
    ext_defaults = read_yaml(find_input_file("extension_defaults.yaml", name = "extension defaults", add_paths=['utils']))
    ignore_units = config.get('ignore_units', False)
    # check_inputs
    all_rules = []
    all_extensions = {'read':{}, 'base':{}, 'other':{}}
    ### get all available inputs instead of just the ones in this config file###
    rules_dir = config['elvers_directories'].get('rules', 'rules')
    all_rules = glob.glob(os.path.join( rules_dir, '*', '*.rule'))
    all_params = {}
    rulenames = []
    for rule in all_rules:
        try:
            rule_name = os.path.basename(rule).split('.rule')[0]
            rulenames.append(rule_name)
            all_params[rule_name] = get_params(rule_name, os.path.dirname(rule))
        except:
            sys.stderr.write(f"\n\tError: Can't decipher params.yml for elvers rule {rule_name}. Please fix.\n\n")
    available_exts = {}
    for rulename, val in all_params.items():
        output_options = val['elvers_params']['output_options']
        #exts = val['elvers_params']['output_options']
        for out_name, val in output_options.items():
            if val['extensions'].get('read'):
                all_extensions['read'][out_name] =  val
            if val['extensions'].get('base'):
                all_extensions['base'][out_name] = val
            if val['extensions'].get('other'):
                all_extensions['other'][out_name] = val

    available_exts['available_extensions'] = all_extensions
    utils_dir = os.path.dirname(os.path.abspath(__file__))
    write_yaml(available_exts, os.path.join(utils_dir, 'extension_defaults.yaml'))
    include_rulenames = [os.path.basename(x).split('.rule')[0] for x in config['include_rules']]
    for rule in rulenames:
        if rule in config.keys():
            if rule not in include_rulenames:
                del config[rule] # bc we won't have built all the elvers_params.
            else:
                config[rule] = generate_rule_targs(home_outdir, base, ref_exts, config[rule], rule, samples, all_extensions, ignore_units)
    return config

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
    paramsD = read_yaml(rule_paramsfile)
    rule_params = paramsD[rule_name]
    return rule_params
