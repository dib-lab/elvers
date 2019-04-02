#! /usr/bin/env python

"""
Execution script for snakemake elvers.
"""
# ref: https://github.com/ctb/2018-snakemake-cli/blob/master/run
import argparse
import os
import sys
import pprint
import yaml
import glob
import collections
import snakemake
import shutil
import subprocess

from snakemake.utils import validate

from .utils.utils import *
from .utils.pretty_config  import pretty_name, write_config
from .utils.print_workflow_options import print_available_workflows_and_tools
from .utils.capture_stdout import CaptureStdout
from .utils.generate_yaml_schema import *

from . import _program

def build_default_params(workdir, targets):
    defaultParams = {}
    # first, figure out which parts of the pipeline are being run, and get those defaults
    pipeline_defaultsFile = find_input_file(os.path.join('utils', 'pipeline_defaults'), 'pipeline_defaults', add_paths = [workdir])
    pipeline_defaults = read_yaml(pipeline_defaultsFile)
    # grab general defaults
    defaultParams['basename'] = pipeline_defaults['basename']
    # add main directories, outdirs available to all workflows
    defaultParams['elvers_directories'] = pipeline_defaults['elvers_directories']
    # grab all available workflows, and subset by user input
    ep_workflows = pipeline_defaults['elvers_workflows']
    workflows_to_run = {k: ep_workflows[k] for k in ep_workflows.keys() & targets}
    defaultParams['elvers_workflows'] = workflows_to_run
    # find targs the user entered that are not in our default info.
    extra_targs = [t for t in targets if t not in ep_workflows]
    for e in extra_targs: # assume extra targets are single rules, and add to the workflow
        workflows_to_run[e] = {'include': [e], 'targets': [e]}
    # For each rule in the desired workflow, save rulename and grab rule params
    required_rules = []
    for targD in workflows_to_run.values():
        required_rules+= targD.get('include', [])
        required_rules+= targD.get('targets', [])
    ruleParamsFiles = []
    includeRules = []
    reference_extensions = []
    rules_dir = defaultParams['elvers_directories']['rules']
    required_rules = set(required_rules)
    for rule_name in required_rules:
        try:
            rule = glob.glob(os.path.join(workdir, rules_dir, '*', rule_name + '.rule'))[0]
            defaultParams[rule_name] = get_params(rule_name, os.path.dirname(rule))
            ref_exts = defaultParams[rule_name]['elvers_params']['outputs']['extensions'].get('reference_extensions', [])
            reference_extensions+=ref_exts
            includeRules += [rule]
        except: # if allows user workflow input, can't do this here (check extra targs later?)
            sys.stderr.write(f"\n\tError: Can't add rules for extra target {rule_name}. Please fix.\n\n")
            sys.exit(-1)
    defaultParams['include_rules'] = list(set(includeRules))
    defaultParams['reference_extensions'] = list(set(reference_extensions))
    return defaultParams


def build_dirs(ep_dir, params):
    ''' function to build full paths for all directories '''
    # build main elvers dir info
    params['elvers_directories']['base_dir'] = ep_dir
    params['elvers_directories']['rules'] = os.path.join(ep_dir, params['elvers_directories']['rules'])
    params['elvers_directories']['animals'] = os.path.join(ep_dir, params['elvers_directories']['animals'])
    # if desired, user can also provide out_path, and all dirs will be built under there
    out_path = params.get('out_path', os.getcwd())
    out_path = os.path.expanduser(out_path) # expand any `~` on unix
    if os.path.isabs(out_path):  # if user inputs an absolute path, check that it exists!
        assert os.path.exists(out_path) and os.path.isdir(out_path), f"Error: provided output path {out_path} is not an existing directory. Please fix.\n\n"
    else: # if not absolute, assume subdirectory of base elvers dir
        out_path = os.path.join(ep_dir, out_path)
    # allow user to define basename, and experiment, build outdir name
    basename = params['basename']
    expt = params.get('experiment', '')
    if expt and not expt.startswith('_'):
        expt= '_' + expt
    outdir = basename + "_out" + expt

    # Now join out_path, outdir name
    out_path = os.path.realpath(out_path)
    outdir = os.path.join(out_path, outdir)
    # when using out_path, need to manually build the outdir (snakemake will not automatically create it)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # add full path info to the config
    params['elvers_directories']['out_dir'] = outdir # outdir NAME
    params['elvers_directories']['logs'] = join(outdir, params['elvers_directories']['logs'])

    # build dirs for main elvers output directories
    outDirs = params['elvers_directories']['outdirs']
    for targ, outD in outDirs.items():
        outDirs[targ] = os.path.join(outdir, outD)
    # put joined paths back in params file
    params['elvers_directories']['outdirs'] = outDirs
    # build dirs for included rules
    included_rules = params['include_rules']
    for rule in included_rules:
        prog = os.path.basename(rule).split('.rule')[0]
        # if no outdir, just use program name
        prog_dir = params[prog]['elvers_params']['outputs'].get('outdir', prog)
        params[prog]['elvers_params']['outputs']['outdir'] = os.path.join(outdir, prog_dir)
    return params

def main():
    parser = argparse.ArgumentParser(prog = _program, description='run snakemake elvers', usage='''elvers <configfile.yaml>  [<target> ...]

Run elvers snakemake workflows, using the given configfile.

Available Workflows:

   default           - run full eel_pond workflow
   preprocess         - preprocess reads
   assemble           - transcriptome assembly
   annotate           - annotate transcriptome assembly
   quantify           - read quantification
   diffexp            - conduct differential expression

For a quickstart, run this:

    elvers examples/nema.yaml

from the main elvers directory.

To build an editable configfile to start work on your own data, run:

    elvers my_config --build_config


''')

    parser.add_argument('configfile')
    parser.add_argument('targets', nargs='*', default=['default'])
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('--extra_config', action='append', default = None)
    parser.add_argument('--config_dict', type=yaml.safe_load, default = None)
    parser.add_argument('--out_path', type=str, default = None)
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-w', '--print_workflows', action='store_true', help='just show available workflows')
    parser.add_argument('-r', '--print_rules', action='store_true', help='just show available rules')
    parser.add_argument('-p', '--print_params', action='store_true', help='print parameters for chosen workflows or rules')
    parser.add_argument('--build_config', action='store_true', help='just build the default parameter file')

    # advanced args below (maybe separate so these don't always print out)
    parser.add_argument('--report', default="report.html", help='filename for a final report of this run. This will be in the logs dir, unless you provide an absolute path.')
    parser.add_argument('--conda_prefix', default=None, help='location for conda environment installs')
    parser.add_argument('--create_envs_only', action='store_true', help="just install software in conda envs, don't execute workflows")
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--dag', action='store_true', help='boolean: output a flowchart of the directed acyclic graph of this workflow in graphviz dot format (does not require/run dot)')
    parser.add_argument('--dagfile', default=None, help='filename to output a flowchart of the directed acyclic graph of this workflow in graphviz dot format')
    parser.add_argument('--dagpng',  default=None, help='filename to output a flowchart of the directed acyclic graph of this workflow in png image format')
    parser.add_argument('-k', '--keep_going', action='store_true')
    parser.add_argument('--nolock', action='store_true')
    parser.add_argument('--unlock', action='store_true')
    parser.add_argument('--cleanup_conda', action='store_true')
    parser.add_argument('--forcetargets', action='store_true', help='force given targets to be re-created (default False)')
    parser.add_argument('--forceall', action='store_true', help='force all output files to be re-created (default False)')
    parser.add_argument('--restart_times',type=int, default=0, help='let snakemake try rerunning any failed tools (input number of times to try rerunning). default = 0')
    args = parser.parse_args()

    thisdir = os.path.abspath(os.path.dirname(__file__))
    # print available workflows and rules, if desired
    if args.print_workflows or args.print_rules:
        pipeline_defaultsFile = find_input_file(os.path.join('utils', 'pipeline_defaults'), 'pipeline_defaults', add_paths = [thisdir])
        print_available_workflows_and_tools(read_yaml(pipeline_defaultsFile), args.print_workflows, args.print_rules)
        sys.exit(0)
    targs = args.targets
    if args.print_params:
        default_params = build_default_params(thisdir, targs)
        write_config(default_params, targs)
        sys.exit(0)

    # are we building a directed acyclic graph?
    building_dag = False
    if args.dag or args.dagfile or args.dagpng:

        building_dag = True

        # building dags precludes running any workflows
        args.dry_run = True

        # if user specified --dagpng,
        # graphviz dot must be present
        if args.dagpng:
            if shutil.which('dot') is None:
                sys.stderr.write(f"\n\tError: Cannot find 'dot' utility, but --dotpng flag was specified. Fix this by installing graphviz dot.\n\n")
                sys.exit(-1)

    # first, find the Snakefile and configfile
    if not building_dag:
        print('\n--------')
        print('checking for required files:')
        print('--------\n')

    snakefile = find_Snakefile(thisdir)

    if args.build_config:
        configfile = args.configfile
        if not any(ext in args.configfile for ext in ['.yaml', '.yml']):
            configfile = args.configfile + '.yaml'
        if os.path.exists(configfile):
            sys.stderr.write(f"\n\tError: found configfile path at {configfile}, but you have '--build_config' specified. Please fix.\n\n")
            sys.exit(-1)
        default_params = build_default_params(thisdir, targs)
        write_config(default_params, targs, configfile)
        sys.exit(0)
    else:
        configfile = find_input_file(args.configfile, 'configfile', add_paths=[thisdir], verbose=True) # find configfile
        if not configfile:
            sys.stderr.write('Error: cannot find configfile {}\n.'.format(args.configfile))
            sys.exit(-1)
        # first, grab all params in user config file
        configD = import_configfile(configfile)
        if configD.get('workflows', None):
            # how do we want to handle the 'default'? Here: If nothing specified, use `default`. If any workflows specified (commandline or config), do not add default.
            if targs == ['default']:
                targs = configD['workflows']
            else:
                targs = targs + configD['workflows']
        # build info for get_reference
        refInput = configD.get('get_reference', None)
        if refInput:
            targs+=['get_reference']
            configD, refinput_ext = handle_reference_input(configD, configfile)
        else:
            refinput_ext = None
        if 'get_reference' in targs and not refInput:
            sys.stderr.write("\n\tError: trying to get reference via `get_reference` rule, but there's no reference file specified in your configfile. Please fix.\n\n")
            sys.exit(-1)
        # check that samples file exists, targs include get_data, and build fullpath to samples file
        samples = None
        if configD.get('get_data', None):
            targs+=['get_data']
            try:
                configD = handle_samples_input(configD, configfile)
                samples = read_samples(configD)
            except Exception as e:
                sys.stderr.write("\n\tError: trying to get input data, but can't find the samples file. Please fix.\n\n")
                print(e)
                sys.exit(-1)
        targs = list(set(targs))
        # next, grab all elvers defaults, including rule-specific default parameters (*_params.yaml files)
        paramsD = build_default_params(thisdir, targs)
       ###############
        # Handle additional configuration modification
        # 1. extra config files
        extra_configs = {}
        if args.extra_config:
            for c in args.extra_config:
                extra_configs = import_configfile(find_input_file(c, 'extra_config', add_paths = [thisdir], verbose =True), extra_configs)

        # 2. config_dict passed in on command line
        # ADVANCED ONLY - no checks in place, formatting matters. (to do: add checks)
        if args.config_dict:
            # adding this the same way as an extra config yaml file
            update_nested_dict(extra_configs, args.config_dict)

        # update_nested_dict only updates keys that already exist. so we need to wait till here (after importing params) to read in extra_configs,
        # rather than doing this when reading in the main configfile

        # if we're adding any additional contrasts, get rid of the example ones!
        # check extra configs
        if extra_configs.get('deseq2'):
           if extra_configs['deseq2']['program_params'].get('contrasts'):
               paramsD['deseq2']['program_params']['contrasts'] = {} # get rid of default contrasts

        # check main configfile
        if configD.get('deseq2'):
            if configD['deseq2']['program_params'].get('contrasts'):
                paramsD['deseq2']['program_params']['contrasts'] = {} # get rid of default contrasts

        # first update with extra configs, then with main configfile
        update_nested_dict(paramsD,extra_configs)

        # update defaults with user-specified parameters (main configfile)
        update_nested_dict(paramsD,configD) # configD takes priority over default params

        # add extension to overall reference_extensions info
        if refinput_ext: # note, need to do it here to prevent override with defaults
            paramsD['reference_extensions'] = list(set(paramsD.get('reference_extensions', []) + [refinput_ext]))

        # This is now handled in the deseq2 rule, can remove from here. Do we need the sys.stderr notification?
        if paramsD.get('no_gene_trans_map', False):
            if paramsD.get('deseq2'):
                paramsD['deseq2']['program_params']['gene_trans_map'] = False
                sys.stderr.write("\tYou're using `get_reference` without specifying a gene-trans-map. Setting differential expression to transcript-level only. See https://dib-lab.github.io/elvers/deseq2/for details.\n")
        # use params to build directory structure
        paramsD = build_dirs(thisdir, paramsD)

        # okay, now lets generate targs and check all required inputs are available
        paramsD = generate_inputs_outputs(paramsD, samples)

        # aggregate the yaml schema for the pipeline (and all included rules) so we can validate against it
        schemafile = os.path.join(os.path.dirname(configfile), '.ep_' + os.path.basename(configfile).rsplit('.y')[0] + '.schema.yaml')
        rulenames = [os.path.basename(x).split('.rule')[0] for x in paramsD['include_rules']]
        pipeline_schema = build_params_schema(paramsD, schemafile, rules=rulenames, targets =targs)
        try:
            # validate the params dictionary we've built using config schema generated for the included targets, all rules paramss
            validate(paramsD, schema=schemafile)
        except Exception as e:
           print(e)
           sys.exit(-1)
        #import pdb;pdb.set_trace()
        #try:
        #    check_workflow(paramsD)
        #except Exception as e:
        #    print(e)
        #    sys.exit(-1)
        # Note: Passing a configfile allows nested yaml/dictionary format.
        # Passing these params in via `config` would require a flat dictionary.
        paramsfile = os.path.join(os.path.dirname(configfile), '.ep_' + os.path.basename(configfile))
        sys.stderr.write('\tAdded default parameters from rule-specific params files.\n\tWriting full params to {}\n'.format(paramsfile))
        write_yaml(paramsD, paramsfile)
        reportfile = None
        if args.report:
            if os.path.isabs(args.report):
                reportfile = args.report
            else:
                reportfile = os.path.abspath(os.path.join(paramsD['elvers_directories']['logs'], args.report))

        if not building_dag:
            print('--------')
            print('details!')
            print('\tsnakefile: {}'.format(snakefile))
            print('\tconfig: {}'.format(configfile))
            print('\tparams: {}'.format(paramsfile))
            print('\ttargets: {}'.format(repr(targs)))
            print('\toutput: {}'.format(repr(paramsD['elvers_directories']['out_dir'])))
            print('\treport: {}'.format(repr(reportfile)))
            print('--------')


        # Set up a context manager to capture stdout if we're building
        # a directed acyclic graph (which prints the graph in dot format
        # to stdout instead of to a file).
        # If we are not bulding a dag, pass all output straight to stdout
        # without capturing any of it.
        passthru = not building_dag
        with CaptureStdout(passthru=passthru) as output:
            # run!!
            # params file becomes snakemake configfile
            status = snakemake.snakemake(snakefile, configfile=paramsfile, use_conda=True,
                                     targets=['elvers'], printshellcmds=True,
                                     cores=args.threads, cleanup_conda= args.cleanup_conda,
                                     dryrun=args.dry_run, lock=not args.nolock,
                                     unlock=args.unlock,
                                     verbose=args.verbose, debug_dag=args.debug,
                                     conda_prefix=args.conda_prefix,
                                     create_envs_only=args.create_envs_only,
                                     restart_times=args.restart_times,
                                     printdag=building_dag, keepgoing=args.keep_going,
                                     forcetargets=args.forcetargets,forceall=args.forceall)

        if building_dag:

            # These three --build args are mutually exclusive,
            # and are checked in order of precedence (hi to low):
            # --dag         to stdout
            # --dagfile     to .dot
            # --dagpng      to .png

            if args.dag:
                # straight to stdout
                print("\n".join(output))

            elif args.dagfile:
                with open(args.dagfile,'w') as f:
                    f.write("\n".join(output))
                print(f"\tPrinted workflow dag to dot file {args.dagfile}\n\n ")

            elif args.dagpng:
                # dump dot output to temporary dot file
                with open('.temp.dot','w') as f:
                    f.write("\n".join(output))
                subprocess.call(['dot','-Tpng','.temp.dot','-o',args.dagpng])
                subprocess.call(['rm','-f','.temp.dot'])
                print(f"\tPrinted workflow dag to png file {args.dagpng}\n\n ")

        if status and reportfile and not args.dry_run and not args.unlock and not building_dag and not args.create_envs_only:
            snakemake.snakemake(snakefile, configfile=paramsfile, report=reportfile)
            print(f"\t see the report file at {reportfile}\n\n ")

        if status: # translate "success" into shell exit code of 0
           return 0
        return 1


# TODO: would be good to pull available rules from elvers_pipeline in default config or a separate workflow file!
if __name__ == '__main__':
    sys.exit(main())
