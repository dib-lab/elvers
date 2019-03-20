import os
import sys
import yaml
import argparse
## generate json schema from a yaml file

# need to disable aliasing for the "properties" dictionaries to dump the appropriate yaml schema
noalias_dumper = yaml.dumper.SafeDumper
noalias_dumper.ignore_aliases = lambda self, data: True

def read_yaml(filename):
    with open(filename, 'r') as stream:
        try:
            yamlD = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return yamlD

def write_yaml(yamlD, paramsfile):
    with open(paramsfile, 'w') as params:
        yaml.dump(yamlD, stream=params, indent=2,  default_flow_style=False, Dumper=noalias_dumper)


def build_program_params_schema(program, params):
    properties = {}
    program_params = params[program]['program_params']
    for key, val in program_params.items():
        if isinstance(val, str):
            properties[key] = {"type": "string"}
        elif isinstance(val, list):
            properties[key] = {"type": "array"}
        elif isinstance(val, bool):
            properties[key] = {"type": "boolean"}
        elif isinstance(val, (int, float)):
            properties[key] = {"type": "number"}
        elif isinstance(val, dict):
            properties[key] = {"type": "object"}
        else:
            print(f"Cannot determine type of {key} default input, {val}")
            pass
    return properties


def build_yaml_schema(args):
    template = read_yaml(args.template)
    params = read_yaml(args.params)
    program_names = list(params.keys())
    final_schema = template.copy()
    final_schema['properties'] = {}
    for name in program_names:
        prog_schema = template['properties']['__prog__'].copy()
        prog_schema['properties']['program_params']['properties'] = build_program_params_schema(name, params)
        final_schema['properties'][name] = prog_schema
    final_schema['required'] = program_names
    write_yaml(final_schema, args.schema_output)


if __name__ == '__main__':
    """ very simple attempt at auto generating yaml schema for elvers rules"""
    p = argparse.ArgumentParser()
    p.add_argument('params')
    p.add_argument('--template', default= "../schemas/param_schema_template.yaml")
    p.add_argument('-o', '--schema_output', default = 'testing_schema.yaml') #default = args.params.split('.')[0] + '_schema.yaml')
    args = p.parse_args()
    sys.exit(build_yaml_schema(args))
