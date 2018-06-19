import os
import yaml
import pandas as pd

def is_se(units,sample, unit, end = ''):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_params(rule_name, rules_dir='rules'):
    rule_paramsfile = os.path.join(rules_dir, rule_name, rule_name + '_params.yaml') # or maybe glob for yaml file in the rule subdir
    rule_params = {}
    with open(rule_paramsfile, 'r') as stream:
        try:
            paramsD = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
        rule_params= paramsD[rule_name]
    return rule_params


