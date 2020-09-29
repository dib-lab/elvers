import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir
import elvers.utils.utils as ep


ep.generate_dir_fullpaths(config)
out_dir = config["output_dir"]
data_dir = config["elvers_directories"]["input_data"]
logs_dir = config["elvers_directories"]["logs"]
benchmarks_dir = config["elvers_directories"]["benchmarks"]
basename = config["basename"]

strict_val = config.get('strict', '1')
strict_mode = int(strict_val)
if not strict_mode:
    print('** WARNING: strict mode is OFF. Config errors will not force exit.')

force = config.get('force', '0')
force = int(force)
force_param = ''
if force:
    force_param = '--force'

# snakemake workflow

# note, this function *needs* to be in this file, or added somewhere it can be accessed by all rules
def is_single_end(name):
    return pd.isnull(samples.at[name, "fq2"])

documentation_base = "https://dib-lab.github.io/elvers/"

octopus = srcdir("utils/animals/octopus")
failwhale = srcdir("utils/animals/failwhale")
fish = srcdir("utils/animals/fish")

onstart:
    print("----------------------------------------------------------")
    print("Elvers, a system for conducting  de novo RNA-Seq Analyses")
    print("----------------------------------------------------------")
    shell('cat {fish}')


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")
    shell('cat {octopus}')
    
    print(f"Outputs for all workflow steps can be found in {out_dir}:\n")
    pipeline = config["pipeline"]
    #print(f"Outputs can be found in {out_dir}")
    steps = config["elvers_pipelines"][pipeline]["steps"]
    for step in steps:
        step_dir = os.path.join(out_dir, config[step]["output_dir"])
        print(f"\t{step}: {step_dir}")
        step_docs = documentation_base + step
        print(f"\t\tdocs: {step_docs} \n\n")

onerror:
    print("  Oh no! Something went wrong here\n")
    shell('cat {failwhale}')

rule elvers:
    input: ep.generate_targets(config)


include_rules = config["include_rules"]
for r in include_rules:
    include: r

# check config files only
rule check:
    input:
        config["sample_info"]                  # do nothing - this file should exist

# print out the configuration
rule showconf:
    input:
        config["sample_info"]
    run:
        import yaml
        import json
        print('# full aggregated configuration:')
        generalP = {}
        generalP['basename'] = config["basename"]
        generalP['experiment_suffix'] = config['experiment_suffix']
        pipeline = config["pipeline"]
        generalP['pipeline'] = pipeline
        step_params = {}

        for step in config["elvers_pipelines"][pipeline]["steps"]:
            # json hack to convert ordered dict to regular dict
            param_dict = json.loads(json.dumps(config[step]["params"]))
            step_params[step] = param_dict
        generalP['program_params'] = step_params
        print(yaml.safe_dump(generalP, indent=2, default_flow_style=False).strip())

        print('# END')
