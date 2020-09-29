import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir
import elvers.utils.utils as ep

out_dir = os.path.abspath(config["output_dir"])
ddir = config["get_data"]["output_dir"]
data_dir = os.path.join("output_dir", ddir)
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")
#report_dir = os.path.join(out_dir, "reports")
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
def is_single_end(sample, end = '', assembly = ''):
    return pd.isnull(samples.at[sample, "fq2"])

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
        print('# full aggregated configuration:')
        generalP = {}
        generalP['basename'] = config["basename"]
        generalP['experiment_suffix'] = config['experiment_suffix']
        print(yaml.dump(generalP, indent=2,  default_flow_style=False).strip())
        pipeline = config["pipeline"]
        print(f"Workflow: {pipeline}")
        step_params, all_program_params = {},{}

        for step in config["elvers_pipelines"][pipeline]["steps"]:
            step_params[step] = config[step]["params"]

        all_program_params["program_params"] = step_params
        print(yaml.dump(all_program_params, indent=2,  default_flow_style=False).strip())
        #print(yaml.dump(all_program_params).strip())
        print('# END')
