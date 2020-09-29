import os
import sys
import yaml
import collections
import pandas as pd
from os.path import join
from snakemake.utils import validate
from snakemake.workflow import srcdir

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
    samples_file = find_input_file(config["sample_info"], "sample info", add_suffixes = ['.tsv'. '.csv', '.xls', '.xlsx'])
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            #samples = pd.read_csv(samples_file, dtype=str, sep=separator).set_index(["sample", "unit"], drop=False)
            samples = pd.read_csv(samples_file, dtype=str, sep=separator).set_index(["sample"], drop=False)
            #validate(samples, schema=srcdir("schemas/samples_v2.schema.yaml"))
            #samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str).set_index(["sample"], drop=False)
            #samples = pd.read_excel(samples_file, dtype=str).set_index(["sample", "unit"], drop=False)
            #validate(samples, schema=srcdir("schemas/samples_v2.schema.yaml")) # new version
            #samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    # column 4 is "condition", but can change name --> LOCK INTO CONDITION
    if (samples.iloc[:, 4].value_counts() < 2).any():
        config['all_replicated'] = False
    # now verify they exist
    data_dir = config.get("data_dir", "")
    if data_dir:
        data_dir = sanitize_path(data_dir)
    #sample_list = samples["fq1"].tolist() + samples["fq2"].tolist()
    for sample in samples.index:
        fq1 = samples.at[sample, "fq1"]
        fullpath_fq1 = os.path.join(data_dir, fq1)
        fq2 = samples.at[sample, "fq2"]
        fullpath_fq2 = os.path.join(data_dir, fq2)
        if not os.path.exists(fullpath_fq1):
            print(f'** ERROR: {sample} fastq1 file {fq1} does not exist in {data_dir}')
            if strict_mode:
                print('** exiting.')
                sys.exit(-1)
            else:
                print(f'** Strict mode is off. Removing this sample and proceeding.')
                samples.drop(sample, inplace=True)
        elif fq2 and not os.path.exists(fullpath_fq2):
            print(f'** ERROR: sample {sample} fastq2 file {fq2} does not exist in {data_dir}')
            if strict_mode:
                print('** exiting.')
                sys.exit(-1)
            else:
                print(f'** Strict mode is off. Treating this sample as single-end and proceeding.')
                samples.at[sample, "fq2"] = ""

    config["samples"] = samples


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

def handle_references(config, pipeline):
    references = []
    gtmaps = []
    input_reference = config.get("input_reference", None)
    urls_begin = config["urls_begin"]
    ref_rule = srcdir("rules/utils/get_reference.rule")
    # handle user-input reference (singular)
    if input_reference and not input_reference.startswith(tuple(urls_begin):
        # check that the input file exists
        input_reference = find_input_file(input_reference, name="input reference", add_suffixes = ['.fa', '.fasta'])
        references.append("refinput")
        config["include_rules"] = [ref_rule]
    gtmap = config.get('input_gene_trans_map', '')
    if gtmap:
        # check that the input file exists
        gtmap = find_input_file(gtmap,"input reference gene_trans_map", add_suffixes = [''])
        gtmaps.append("refinput")

    # handle assembler-generated references
    steps = config["elvers_workflows"][pipeline]["steps"]
    for step in steps:
        if config[step].get("generates_reference"):
            references.append(step)
    # make sure at least one ref or assembly program exists!
    if not references:
        sys.stderr.write("\n\t Please input a transcriptome reference via 'input_reference` or include an assembly program (e.g. trinity, plass) in the workflow.\n\n")
        sys.exit(-1)
    refdir = config["get_reference"]["output_dir"]
    reference_targets = config["get_reference"]["output_files"]
    genetransmap_targets = config["get_reference"]["output_gene_trans_map"]

    #add reference names to config
    config["references"] = references
    output_dir= config["output_dir"]
    # build reference targets
    ref_targets = expand(os.path.join(output_dir, refdir, reference_targets, basename = config["basename"], reference = references)
    gtm_targets = expand(os.path.join(output_dir, refdir, genetransmap_targets, basename = config["basename"], reference = references)

    return ref_targets + gtm_targets

# sample checks
def is_single_end(sample, end = '', assembly = ''):
    return pd.isnull(samples.at[sample, "fq2"])

def handle_samples_input(config):
    program_params = config['get_data'].get('program_params')
    urls_begin = config["urls_begin"]
    data_rule = srcdir("rules/utils/get_data.rule")
    samples_file = config.get('sample_info', None)
    if samples_file:
        read_samples(samples_file)
        config["include_rules"].append(data_rule)
    else:
        sys.stderr.write("\n\tError: this workflow needs samples files, but no 'samples_info' file is provided in the configfile. Please fix.\n\n")
        sys.exit(-1)

def handle_user_program_params(config):
    # snakemake doesn't properly handle nested configuration dictionaries
    # handle manually instead:
    user_program_params = config.get("program_params", {})
    for program, user_params in user_program_params.items():
        if program not in config.keys():
            sys.stderr.write(f"\nWarning: New parameters for program {program} provided in the configfile, but this program name doesn't match any in this pipeline. Ignoring.\n\n")
             continue
        params = config[program].get("params", {})
        update_nested_dict(params, user_params)
        config[program]["params"] = params


# make sure we have sample and /or reference info
def generate_targets(config):
    # get info for pipeline we're running
    pipeline = config["pipeline"]
    basename = config["basename"]
    reference_targets, workflow_targets=[],[]
    config["include_rules"] = []

    if config["elvers_pipelines"][pipeline]["reference_required"]:
        reference_targets = handle_references(config)
    if config["elvers_pipelines"][pipeline]["samples_required"]:
        handle_samples_input(config)
        workflow_targets = generate_pipeline_targets(config, pipeline, config["samples"])

    return reference_targets + workflow_targets

def generate_pipeline_targets(config, pipeline, samples):
    pipeline_targets=[]
    pipeline_rules = []
    # generate targets for each step
    steps = config["elvers_pipelines"][pipeline]["steps"]
    config["workflow_steps"] = steps

    if samples:

    for step in steps:
        rulefile = config[step]["rulefile"]
        r = srcdir(f"rules/{rulefile}")
        config['include_rules'].append(r)

        step_outdir = config[step]["output_dir"]
        step_files = config[step]["output_files"]

        references = config["references"]
        pe_ends = config['fq_ends']["paired"]
        pe_pairing = config['pairing']["paired"]
        contrasts = ""

        for stepF in step_files:
            if stepF["params"].get("contrasts", None):
                contrasts = stepF["params"]["contrasts"].keys()
            for sample in samples:
                # "ends", ""pairing" will be different for se, pe files
                pairing = pe_pairing
                ends = pe_ends
                if is_single_end(sample, end = '', assembly = ''):
                    pairing = se_pairing
                    ends = se_ends
                pipeline_targets += expand(os.path.join(output_dir, step_outdir, stepF), sample=sample, reference=references, basename=basename, pairing=pairing, end=ends, contrast=contrasts)

    return pipeline_targets


