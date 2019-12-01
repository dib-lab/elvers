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

def read_samples(config, build_sra_links = False):
    elvers_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # find samples file
    try:
        samples_file = config["get_data"]["program_params"]["samples"]
    except KeyError:
        print("cannot find 'samples' entry in config file! Please provide a samples file for the `get_data` utility", file=sys.stderr)
        sys.exit(-1)
    # read samples file
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            samples = pd.read_csv(samples_file, dtype=str, sep=separator).set_index(["sample", "unit"], drop=False)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str).set_index(["sample", "unit"], drop=False)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    # build sra links if necessary (reads)
    if build_sra_links:
        base_link = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
        if "SRR" in samples.columns and "LibraryLayout" in samples.columns:
            samples['fq1'] = samples['SRR'].apply(lambda x : base_link + x[0:6] + '/00' + x[-1] + '/' + x + '/' + x + '_1.fastq.gz')
            samples['fq2'] = df.apply(lambda row : build_fq2(row), axis=1)
    # check ("validate") that the required columns are present
    try:
        validate(samples, schema=os.path.join(elvers_dir,"schemas/samples_v2.schema.yaml"))
    except Exception as e:
        sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
        print(e)
    samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
    # check for single-unit case
    if not (samples['sample'].value_counts() > 1).any():
        config['ignore_units'] = True
    ## now, we restrict replicate checking to "condition" name. Make a note in docs
    if 'condition' in samples.columns:
        if (samples['condition'].value_counts() < 2).any():
            config['all_replicated'] = False
    else:
        config['all_replicated'] = False
    return samples, config

def build_fq2(row):
    base_link = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    SRR = row['SRR']
    if row['LibraryLayout'] == "PAIRED":
        fq2 = base_link + SRR[0:6] + '/00' + SRR[-1] + '/' + SRR  + '/' + SRR + '_2.fastq.gz'
    return fq2

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

def handle_reference_input(config, configfile, samples = None):
    reference_extensions = []
    program_params = config['get_reference'].get('program_params')
    firstref = {}
    reference_list = {}

    # current, single ref specification allows no reference extension
    if program_params.get('reference', None):
        initial_ref = {'reference': program_params['reference'],   \
                       'gene_trans_map': program_params.get('gene_trans_map', None), \
                       'reference_extension': program_params.get('reference_extension', ""), \
                       'associated_samples': program_params.get('associated_samples', None)}

        program_params = check_ref_input(initial_ref, configfile) #check for file inputs & if yes, return with fullpaths
        ref_ext = initial_ref['reference_extension'] #program_params.get('reference_extension', "")
        if ref_ext:
            if not ref_ext.startswith('_'):
                ref_ext = '_' + ref_ext
            reference_extensions = [ref_ext]
            firstref[ref_ext] = initial_ref
        # just trying this out -- would mean we only have to deal with a single ref specification
        else:
            firstref["no_extension"] = initial_ref
            reference_extensions = ["no_extension"]

    # multiple ref specification (each ref needs an extension)
    if program_params.get('reference_list', None):
        reference_list = program_params['reference_list']
        # should look like: {ref_ext: {reference: file/link, associated_samples: [sample_list], gene_trans_map: gtmap file/link}

        for ref_ext, ref_info in reference_list.items():
            #assoc_samples = ref_info.get('associated_samples', None) # should be samples list --> don't need to change this, leave as-is
            if not ref_ext.startswith('_'):
                ref_ext = '_' + ref_ext
            reference_extensions.append(ref_ext)
            ref_info = check_ref_input(ref_info, configfile)

    per_sample_ref = program_params.get('per_sample_reference_files', False)
    #if 'reference' in samples.columns and not per_sample_ref:
    #    sys.stderr.write("\n\t I see a 'reference' column in your samples file, \
         # do you want to use these as per-sample references? If yes, please add \
        # 'per_sample_reference_files: True' to your 'get_reference' directive in your yaml configuration file\n\n")
    if per_sample_ref:
        if 'reference' not in samples.columns:
            sys.stderr.write("\n\t You have 'per_sample_reference_files' set in your config file, but no 'reference' column in your samples file. Please fix.\n\n")
            sys.exit(-1)
        for row in samples.itertuples(index=False):
            if 'gene_trans_map' in samples.columns:
                sample_ref = {'reference': row.reference, 'gene_trans_map': row.gene_trans_map, 'associated_samples': [row.sample]}
            else:
                sample_ref = {'reference': row.reference, 'gene_trans_map': None, 'associated_samples': [row.sample]}
            sample_ref =  check_ref_input(sample_ref, configfile)
            if not row.sample.startswith('_'):
                reference_extensions.append('_' + row.sample)
            else:
                reference_extensions.append(row.sample)
            # REFERENCE_LIST
            reference_list[row.sample] = sample_ref
    extensions = {'base': ['.fasta'], 'reference_extensions': reference_extensions}
    reference_list.update(firstref)
    #program_params['reference_list'] = reference_list
    config['reference_info'] = reference_list
    config['get_reference'] = {'program_params': program_params, 'elvers_params': {'outputs': {'extensions':extensions}}}
    config['reference_extensions'] = reference_extensions
    ## IF WE HAVE NO REFERENCE INFO BY NOW, SOMETHING IS WRONG - WARN USER
    #if not reffile:
        #    sys.stderr.write("\n\tError: improper reference specification in `get_reference`. Please fix.\n\n")
    return config, reference_extensions


def handle_assemblies(config, samples = None, assembly_programs = ['trinity', 'plass']):
    ref_exts = []
    assembly_info = {}
    # handle per-sample assemblies: use _sample_assembler as the reference extension
    per_sample_assemb = config.get('per_sample_assemblies', [])
    for asmb_prog in per_sample_assemb:
        if not config.get(asmb_prog, None):
            sys.stderr.write(f'Error: per_sample_assemblies specifies the {asmb_prog} assembly program for assemblies, but the program is not included in the workflow config. Did you specify a valid assembler?')
            sys.exit(-1)
        for row in samples.itertuples(index=False):
            if not row.sample.startswith('_'):
                ref_ext = '_' + row.sample + '_' + asmb_prog
            else:
                ref_ext = row.sample + '_' + asmb_prog
            assembly_info[ref_ext] = {'associated_samples': [row.sample]}
            ref_exts.append(ref_ext)
    ref_exts = list(set(ref_exts))
    # process the assembly list from main config file
    assembly_list = config.get('assembly_list', None)
    if assembly_list:
        for assemb_name, info in assembly_list.items():
            assemblers = info.get('assemblers', None)
            if 'associated_samples' in info.keys():
                assoc_samples = info.get('associated_samples')
            else:
                assoc_samples = None
            for assembler in assemblers:
                ref_ext = '_'.join(assemb_name, assembler)
                if not ref_ext.startswith('_'):
                    ref_ext = '_' + ref_ext
                if not config.get(assembler, None):
                    sys.stderr.write(f'Error: The config assembly list specifies the {assembler} assembly program for the {ref_ext} assembly, but the program is not included in the workflow config. Did you specify a valid assembler?')
                    sys.exit(-1)
                if ref_ext in assembly_info.keys():
                    sys.stderr.write(f'Error. The reference extension {ref_ext} in being used for more than one assembly, cannot process this entry in the assembly_list config info')
                    sys.exit(-1)
                elif ref_ext in reference_extensions:
                    sys.stderr.write(f'Error. The reference extension {ref_ext} is being used for both a reference and an assembly, cannot process this entry in the assembly_list config info')
                    sys.exit(-1)
                else:
                    if not assoc_samples:
                        assembly_info[ref_ext] = {'assembler': assembler}
                    else:
                        assembly_info[ref_ext] = {'assembler': assembler, 'associated_samples': assoc_samples}
                ref_exts.append(ref_ext)
    # now, grab assembly info from assembler configs
    for assembler in assembly_programs:
        prog_info = config.get(assembler, None)
        if prog_info:
            prog_params = prog_info.get('program_params')
            try:
                ref_ext = prog_info['elvers_params']['output_options']['extensions']['reference_extensions']
            except:
                ref_ext = '_' + assembler
            if 'associated_samples' in prog_params.keys():
                assoc_samples = prog_params.get('associated_samples', None)
            else:
                assoc_samples = None
            # how to handle duplicate, single assembly. Aka user just wants ONE trinity assembly, but specifies it within assembly_list and now we're checking the trinity info.
            if ref_ext in assembly_info.keys():
                prev_info = assembly_info[ref_ext]
                prev_samples = None
                if prev_info.get('associated_samples', None):
                    prev_samples = set(prev_info['associated_samples'])
                if assoc_samples:
                    assoc_samples = set(assoc_samples)
                if not prev_samples == assoc_samples:
                    p_s = '\n\t' + '\n'.join(prev_samples)
                    a_s = '\n\t' + '\n'.join(assoc_samples)
                    sys.stderr.write(f'Error. The reference extension {ref_ext} is already in use for the samples {p_s}, cannot add associated samples {a_s} from the assembler {assembler} config info')
                    sys.exit(-1)
            else:
                if not assoc_samples:
                    # don't want to automatically assume we want a full assembly, if other assemblies are being generated.
                    all_assemblies = assembly_info.keys()
                    if not any (ref_ext in x for x in all_assemblies):
                        assembly_info[ref_ext] = {'assembler': assembler}
                        ref_exts.append(ref_ext)
                else:
                    assembly_info[ref_ext] = {'assembler': assembler, 'associated_samples': assoc_samples}
                    ref_exts.append(ref_ext)

    reference_extensions = config.get('reference_extensions', [])
    ref_exts = list(set(reference_extensions + ref_exts))
    if assembly_info:
        config['assembly_info'] = assembly_info
    if ref_exts:
        config['reference_extensions'] = ref_exts
    return config


def check_ref_input(refDict, configfile):
    # refDict contains at least "reference", and optionally "gene_trans_map" and/or "reference_extension"

    # if file, find ref
    ref = refDict['reference']
    if not ref.startswith('http') and not ref.startswith('ftp'):
        refDict['reference'] = find_input_file(ref, name="input reference", add_paths = [os.path.realpath(os.path.dirname(configfile))], add_suffixes = ['.fa', '.fasta'])

    # if file, find gene_trans_map
    gtmap = refDict.get("gene_trans_map", None)
    if gtmap:
        if not gtmap.startswith('http') and not gtmap.startswith('ftp'):
            refDict['gene_trans_map'] = find_input_file(gtmap,"input reference gene_trans_map", add_paths = [os.path.realpath(os.path.dirname(configfile))], add_suffixes = [''])
    else:
        del refDict['gene_trans_map']

    # check that we actually have associated samples
    assoc_samples = refDict.get("associated_samples", None)
    if not assoc_samples:
        del refDict['associated_samples']
    return refDict
    # extensions are going to need to change based on each file. Handle get_reference as a special case in building targets, don't set the extensions here
        #extensions = {'base': ['.fasta']}

    ### THIS MEANS WE NEED TO CHANGE DESEQ2 GENE-TRANS-MAP  FINDING/IDENTIFICATION
    #config['no_gene_trans_map']= True

def handle_samples_input(config, configfile):
    program_params = config['get_data'].get('program_params')
    samples_file = program_params.get('samples', None)
    if samples_file:
        samples_file = find_input_file(samples_file, "samples tsv", add_paths = [os.path.realpath(os.path.dirname(configfile))], add_suffixes = ['.tsv'])
        config['get_data']['program_params']['samples'] = samples_file
    else:
        sys.stderr.write("\n\tError: trying to run `get_data` workflow, but the samples tsv file is not specified in your configfile. Please fix.\n\n")
    return config

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
                            break
                            # reference_extensions+=ref_exts
                            ## HERE IS WHERE WE ADD REFERENCE EXTENSIONS FOR ASSEMBLIES
                            #--> not doing this anymore. it will be done in the "handle_assemblies" section
                            #for ref_ext in ref_exts:
                            #    if ref_ext not in reference_extensions:
                            #        reference_extensions.append(ref_ext) # first, append

                        ### MAYBE THIS NEEDS TO NOT BE AN ELSE?
                        elif not val['elvers_params'].get('outputs', None):
                            sys.stderr.write(f"Error: cannot find corresponding outputs for inputs {inputs} for rule {key}")
                            sys.exit(-1)
                    val['elvers_params']['outputs'] = outputs
                else:
                    val['elvers_params']['outputs'] = val['elvers_params']['output_options']
                del val['elvers_params']['output_options']
                config[key] = val
            # we also get rid of the output_options, to minimize clutter
    #config['reference_extensions'] = reference_extensions
    return config


def generate_targs(outdir, basename, samples, ref_exts=[''], base_exts = None, read_exts = None, other_exts = None, contrasts = [], ignore_units=False, ref_info = None, assemb_info = None):
    base_targets, read_targs, other_targs = [],[],[]
    ref_pe_ext, ref_se_ext = [],[]
    if read_exts:
        combine_units = read_exts.get('combine_units', False)
        if combine_units and not ignore_units:
            se_names = samples.loc[samples["fq2"].isnull(), 'sample'].tolist()
            pe_names = samples.loc[samples["fq2"].notnull(), 'sample'].tolist()
        else:
            se_names = samples.loc[samples["fq2"].isnull(), 'name'].tolist()
            pe_names = samples.loc[samples["fq2"].notnull(), 'name'].tolist()
        # grab extensions
        pe_ext = read_exts.get('pe', None)
        ref_pe_ext = [x for x in pe_ext if '__reference__' in x ]
        read_only_pe_ext = [x for x in pe_ext if '__reference__' not in x ]
        # se
        se_ext = read_exts.get('se', None)
        ref_se_ext = [x for x in se_ext if '__reference__' in x ]
        read_only_se_ext = [x for x in se_ext if '__reference__' not in x ]
        # build read targets (if don't use reference info)
        if read_only_se_ext and len(se_names)>0:
            read_targs+=[join(outdir, name + e) for e in read_only_se_ext for name in se_names]
        if read_only_pe_ext and len(pe_names) > 0:
            read_targs+=[join(outdir, name + e) for e in read_only_pe_ext for name in pe_names]
    # handle base targets, e.g. refname_ext.fasta or refname_ext.stats
    # build base targets (using reference extensions)
    #if base_exts:
    #if ref_exts:
        #references = [basename + e for e in ref_exts]
    #else:
    #    references = [basename]
    #for refname in references:
    if ref_exts: # if we have any references at all
        for refx in ref_exts:
            thisref_read_targs = []
            if "no_extension" in refx:
                refname = basename
            else:
                refname = basename + refx
            #references.append(refname)
            if base_exts:
                if assemb_info:
                    # handle additional assemblies
                    rxs = [basename + k for k in assemb_info.keys() if refx in k]
                    base_targets += [join(outdir, rn + e) for e in base_exts for rn in rxs]
                else:
                    # single assemblies
                    base_targets += [join(outdir, refname + e) for e in base_exts]
            # if the read targets use reference_info
            if ref_pe_ext or ref_se_ext:
                # HERE, handle associated samples
                if refx.startswith('_'):
                    refx = refx.split('_')[1]
                ## HERE WE NOW HAVE AN ISSUE FOR ASSEMBLERS. --> build ref_list for assemblies first?
                if refx in ref_info.keys():
                    assoc_samples = ref_info[refx].get('associated_samples', None)

                elif refx in assemb_info.keys():
                    assoc_samples = assemb_info[refx].get('associated_samples', None)
                else:
                    assoc_samples = None

                if assoc_samples:
                    # use subset of samples DF to generate sample names
                    assoc_subset = samples[samples['sample'].isin(assoc_samples)]
                    if combine_units and not ignore_units:
                        se_assoc = assoc_subset.loc[assoc_subset["fq2"].isnull(), 'sample'].tolist()
                        pe_assoc = assoc_subset.loc[assoc_subset["fq2"].notnull(), 'sample'].tolist()
                    else:
                        se_assoc = assoc_subset.loc[assoc_subset["fq2"].isnull(), 'name'].tolist()
                        pe_assoc = assoc_subset.loc[assoc_subset["fq2"].notnull(), 'name'].tolist()
                # handle read targets that contain reference info
                    if ref_se_ext and len(se_assoc)>0:
                        thisref_read_targs+=[join(outdir, name + e) for e in ref_se_ext for name in se_assoc]
                    if pe_ext and len(pe_assoc) > 0:
                        thisref_read_targs+=[join(outdir, name + e) for e in ref_pe_ext for name in pe_assoc]
                else: # if no associated samples, assume ALL, and add all as targs.
                    if ref_se_ext and len(se_names)>0:
                        thisref_read_targs+=[join(outdir, name + e) for e in ref_se_ext for name in se_names]
                    if pe_ext and len(pe_names) > 0:
                        thisref_read_targs+=[join(outdir, name + e) for e in ref_pe_ext for name in pe_names]
                # now replace all instances of "__reference__" for this ref
                thisref_read_targs = [t.replace("__reference__", refname) for t in thisref_read_targs]
                # now add these targs to the overall read targets
                read_targs+= thisref_read_targs #[t.replace("__reference__", refname) for t in read_targets] #should return read_targets if nothing to replace

    #handle contrasts within the base targets
    base_targs = []
    if contrasts and base_targets:
        for c in contrasts:
            base_targs += [t.replace("__contrast__", c) for t in base_targets]
    else:
        base_targs = base_targets
    # handle outputs with no name into (e.g. multiqc)
    if other_exts:
        # no extension, just single filename that we don't want to generate for multiple reference extensions
        other_targs = [join(outdir, e) for e in other_exts]
    return base_targs + read_targs + other_targs

# not in use anymore
#def generate_program_targs(configD, samples, basename,ref_exts, contrasts):
    # given configD from each program, build targets
#    outdir = configD['outdir']
#    exts = configD['extensions']
#    if exts.get('reference_extensions'): # this program is an assembler or only works with specific assemblies
#        ref_exts = exts.get('reference_extensions', ['']) # override generals with rule-specific reference extensions
#    targets = generate_targs(outdir, basename, samples, ref_exts, exts.get('base', None),exts.get('read'), exts.get('other'), contrasts)
#    return targets

def generate_rule_targs(home_outdir, basename, ref_exts, rule_config, rulename, samples, all_input_options, ignore_units, reference_info, assembly_info):
    contrasts = rule_config['program_params'].get('contrasts', [])

    # handle input options!
    if rulename == 'get_data':
        samples_file = rule_config['program_params']['samples'] # should be present (validated) prior to here
        rule_config['elvers_params']['input_options'] = {'get_data': {'indir': os.path.dirname(samples_file), 'input_files': [samples_file]}}

    elif rulename == 'get_reference':
        ## Jointly do inputs/outputs for get_reference
        outdir = rule_config['elvers_params']['outputs']['fasta']['outdir']
        #outdir = output_info['outdir']
        ref_output_files = []
        ref_input_files = []

        # multiple reference input
        # NOW WE SHOULD ALWAYS HAVE A REFERENCE_LIST
        ## NEED TO GENERATE IT FOR ASSEMBLIES!? --> ASSEMBLY_INFO

        #if rule_config['program_params'].get('reference_list'):
        for ref_ext, ref_info in reference_info.items(): #rule_config['program_params']['reference_list'].items():
            thisref_exts = ['.fasta']
            ref_input_files.append(ref_info['reference'])
            if ref_info.get('gene_trans_map'):
                ref_input_files.append(ref_info['gene_trans_map'])
                thisref_exts.append('.fasta.gene_trans_map')
            if not ref_ext.startswith('_'):
                ref_ext = '_' + ref_ext
            thisref = [ref_ext]
            ref_output_files += generate_targs(outdir, basename, samples, ref_exts= thisref, base_exts= thisref_exts, ref_info = reference_info, assemb_info = assembly_info)

        rule_config['elvers_params']['input_options'] = {'get_ref': {'indir': "", 'input_files': ref_input_files}}
        rule_config['elvers_params']['outputs']['output_files'] = ref_output_files
        rule_config['elvers_params']['outputs']['outdir'] = outdir
        return rule_config

    else:
        # handle inputs for all other rules
        all_input_exts = {}
        input_options = rule_config['elvers_params'].get('input_options', None)# read, reference, other
        not_found = []
        input_files = None
        for input_type, options in input_options.items():
            for option in options:
                try:
                    info = all_input_options[input_type].get(option)
                    all_input_exts[option] = info
                    indir = os.path.join(home_outdir, info.get('outdir', 'input_data'))
                    in_exts = info['extensions']
                    if in_exts.get('reference_extensions'): # this program is an assembler or only works with specific assemblies
                        ### REFERENCE EXTENSIONS
                        # NEED TO MODIFY FOR MULTIPLE ASSEMBLIES --> contains rather than =?don't replace but subset the reference_exts?
                        in_ref_exts = in_exts.get('reference_extensions', ['']) # override generals with rule-specific reference extensions
                    else:
                        in_ref_exts = ref_exts
                    input_files = generate_targs(indir, basename, samples, ref_exts =in_ref_exts, base_exts=in_exts.get('base', None),read_exts =in_exts.get('read'), other_exts=in_exts.get('other'), ignore_units=ignore_units, contrasts = contrasts, ref_info = reference_info, assemb_info = assembly_info)
                    all_input_exts[option]['input_files'] = input_files #generate_targs(indir, basename, samples, in_ref_exts, in_exts.get('base', None),in_exts.get('read'), in_exts.get('other'), contrasts)
                    all_input_exts[option]['indir'] = indir
                except:
                    not_found.append(option)
            ## THIS ISN'T QUITE WORKING YET
            if len(not_found) == len(options):
                option_list = "\n  " + "\n  ".join(options)
                sys.stderr.write(f"cannot find input files for {rulename}. Please add a target that produces any of the following: {option_list}")
                sys.exit(-1)
        rule_config['elvers_params']['input_options'] = all_input_exts

    # now handle outputs (sometimes multiple outputs)
    output_files = []
    outdir = ""

    for outname, output_info in rule_config['elvers_params']['outputs'].items():
        outdir = output_info['outdir']
        out_exts = output_info['extensions']
        # NEED TO HANDLE THIS DIFFERENTLY FOR PER-SAMPLE REF ASSEMBLY
        if out_exts.get('reference_extensions'): # this program is an assembler or only works with specific assemblies
            out_ref_exts = out_exts.get('reference_extensions', ['']) # override generals with rule-specific reference extensions
        else:
            out_ref_exts = ref_exts
        outputs = generate_targs(outdir, basename, samples, out_ref_exts, out_exts.get('base', None),out_exts.get('read'), out_exts.get('other'), ignore_units=ignore_units, contrasts=contrasts, ref_info=reference_info, assemb_info = assembly_info)
        output_files += outputs
    rule_config['elvers_params']['outputs']['output_files'] = output_files
    rule_config['elvers_params']['outputs']['outdir'] = outdir
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
    # grab info for all references
    reference_info = {}
    if config.get('get_reference'):
        if config.get('reference_info'):
            reference_info = config['reference_info'] #config['get_reference']['program_params']['reference_list']
        else:
            sys.stderr.write(f"\nError: Trying to run the get_reference utility without specifying reference information. Please fix. \n\n")
            sys.exit(-1)
        #reference_info = config['get_reference']['program_params']['reference_list']
    assembly_info = config.get('assembly_info', {})

     # if we have reference extensions, add them to the reference_info dictionary.
    #if ref_exts:
    #    for ref_ext in ref_exts:
    #        if ref_ext not in assembly_info.keys():
    #            reference_info[ref_ext] = {}

    for rule in rulenames:
        if rule in config.keys():
            if rule not in include_rulenames:
                del config[rule] # bc we won't have built all the elvers_params.
            else:
                config[rule] = generate_rule_targs(home_outdir, base, ref_exts, config[rule], rule, samples, all_extensions, ignore_units, reference_info, assembly_info)
    return config

def get_params(rule_name, rule_dir='rules'):
    # pass in a rule name & the directory that contains its paramsfile.
    # Return paramsD
    rule_paramsfile = os.path.join(rule_dir, 'params.yml')
    paramsD = read_yaml(rule_paramsfile)
    rule_params = paramsD[rule_name]
    return rule_params
