#! /usr/bin/env python
"""
Print all available workflow options
"""
import os
import sys
import glob
import yaml
import argparse

from .utils import read_yaml
from .pretty_config import pretty_name

def print_available_workflows_and_tools(paramsD, print_workflows=True, print_rules=False, only_rules=False):
    if only_rules:
        print_workflows = False
    if print_workflows:
        # print available workflows
        sys.stdout.write(pretty_name("Available elvers Workflows") +'\n')
        ep_flows = paramsD.get('elvers_workflows')
        for workflow, tools in ep_flows.items():
            sys.stdout.write('\n  ' + workflow + ':\n\t')
            for t in tools.values():
                sys.stdout.write('\n\t'.join(t) + '\n\t')

    if print_rules:
        # print available rules
        sys.stdout.write(pretty_name("Advanced usage: all available rules") + '\n\t')
        rules_dir = paramsD['elvers_directories'].get('rules', 'rules')
        rules = glob.glob(os.path.join( rules_dir, '*', '*.rule'))
        rule_names = [os.path.basename(r).split('.rule')[0] for r in rules]
        sys.stdout.write('\n\t'.join(rule_names) + '\n\n\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('paramsfile')
    p.add_argument('-w', '--print_workflows', action='store_true', default=True)
    p.add_argument('-r', '--print_rules', action='store_true', default=False)
    p.add_argument('--only_rules', action='store_true', default=False)
    args = p.parse_args()
    # read params here, to enable passing in paramsD instead of paramsfile from run_elvers
    params = read_yaml(args.paramsfile)
    sys.exit(print_available_workflows_and_tools(params, args.print_workflows, args.print_rules, args.only_rules))
