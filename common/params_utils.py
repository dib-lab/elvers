import yaml

def get_config(RULES_DIR, rule_name):
    rule_paramsfile = join(RULES_DIR, rule_name, rule_name + '_params.yaml') # or maybe glob for yaml file in the rule subdir
    rule_params = {}
    with open(rule_paramsfile, 'r') as stream:
        try:
            config_dir = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return rule_params
