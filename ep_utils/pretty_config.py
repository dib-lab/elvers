#! /usr/bin/env python
"""
Print a readable config yaml for user to edit
"""
import os.path
import argparse
import sys
import yaml

def read_yaml(filename):
    with open(filename, 'r') as stream:
        try:
            yamlD = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return yamlD

def pretty_name(targ):
    split = "  ####################  " 
    name = '\n' + split + targ + split + '\n'
    return name

def write_config(paramsD, targets, out):
    
    print('\n\tprinting editable configfile to {}'.format(out))
    with open(out, 'w') as outConfig:
        # write general eelpond pipeline params
        generalP = {}
        generalP['basename'] = paramsD.get('basename', 'eelpond')
        generalP['experiment'] = paramsD.get('experiment', '_experiment1')
        generalP['samples'] = paramsD.get('samples', 'samples.tsv')
        
        outConfig.write((pretty_name("Eelpond Pipeline Configfile")))
        yaml.dump(generalP, stream=outConfig, indent=2,  default_flow_style=False)
        
        # write program-specific parameters
        seen_rules = []
        for targ in targets:
            # grab all rules, their params for target pipeliness
            include_rules = paramsD['eelpond_workflows'][targ].get('include', [])
            targets = paramsD['eelpond_workflows'][targ].get('targets', [])
            all_rules = include_rules+targets
            targ_params = {}
            for r in all_rules:
                if r not in seen_rules:
                    seen_rules+=[r]
                    if paramsD[r].get('program_params', None):
                        targ_params[r] = paramsD[r]['program_params']
            # write to file
            outConfig.write((pretty_name(targ)))
            yaml.dump(targ_params, stream=outConfig, indent=2,  default_flow_style=False)
    print('\n\tdone! Now edit parameters in the {}, and rerun run_eelpond without the "--build_config" option.\n\n'.format(out))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('paramsfile')
    p.add_argument('-o', '--output')
    p.add_argument('-t', '--targets', nargs='+')
    args = p.parse_args()
    assert args.output, "must specify location for output configfile using '-o'"
    params = read_yaml(args.paramsfile)
    sys.exit(write_config(params, args.targets, args.output))
