#! /usr/bin/env python
"""
Print a readable config yaml for user to edit
"""
import os.path
import argparse
import sys
import yaml

from .utils import read_yaml

def pretty_name(targ):
    print_name = '\n#\n# '+ targ + '\n#\n\n'
    return print_name

def write_config(paramsD, targets, out = None):
    if out:
        print('\n\tprinting editable configfile to {}'.format(out))
        outConfig = open(out, 'w')
        # write general elvers pipeline params
        generalP = {}
        generalP['basename'] = paramsD.get('basename', 'elvers')
        generalP['experiment'] = paramsD.get('experiment', '_experiment1')
        outConfig.write((pretty_name("elvers pipeline configuration")))
        yaml.dump(generalP, stream=outConfig, indent=2,  default_flow_style=False)
    else:
        sys.stdout.write('\n\n')

    # write program-specific parameters
    seen_rules = []
    for targ in targets:
       # grab all rules, their params for target pipeliness
        include_rules = paramsD['elvers_workflows'].get(targ, [])
        all_rules = include_rules #+targets
        targ_params = {}
        for r in all_rules:
            if r not in seen_rules:
                seen_rules+=[r]
                if paramsD[r].get('program_params', None):
                   targ_params[r] = paramsD[r]['program_params']
        # write to file
        if out:
            outConfig.write(pretty_name(targ + ' workflow'))
            yaml.dump(targ_params, stream=outConfig, indent=2,  default_flow_style=False)
        else:
            sys.stdout.write(pretty_name(targ + ' workflow'))
            yaml.dump(targ_params, stream=sys.stdout, indent=2,  default_flow_style=False)
    if out:
        print('\n\tdone! Now edit parameters in the {}, and rerun run_elvers without the "--build_config" option.\n\n'.format(out))
        outConfig.close()
    else:
        sys.stdout.write('\n\n')
        #sys.stdout.write('  #######################################################\n\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('paramsfile')
    p.add_argument('-o', '--output')
    p.add_argument('-t', '--targets', nargs='+')
    args = p.parse_args()
    #assert args.output, "must specify location for output configfile using '-o'"
    params = read_yaml(args.paramsfile)
    sys.exit(write_config(params, args.targets, args.output))
