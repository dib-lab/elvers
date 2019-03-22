import os
import sys
import yaml
import glob
import argparse
## generate json/yaml schema from a yaml file
from elvers.utils.utils import update_nested_dict, find_input_file, read_yaml

elvers_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# this disables yaml aliasing, which allows the "properties" dictionaries to dump properly
noalias_dumper = yaml.dumper.SafeDumper
noalias_dumper.ignore_aliases = lambda self, data: True

def write_yaml(yamlD, paramsfile):
    with open(paramsfile, 'w') as params:
        yaml.dump(yamlD, stream=params, indent=2,  default_flow_style=False, Dumper=noalias_dumper)

def build_schema(param_dict):
    properties = {}
    for key, val in param_dict.items():
        if isinstance(val, str):
            properties[key] = {"type": "string"}
        elif isinstance(val, list):
            properties[key] = {"type": "array"}
        elif isinstance(val, bool):
            properties[key] = {"type": "boolean"}
        elif isinstance(val, (int, float)):
            properties[key] = {"type": "number"}
        elif isinstance(val, dict):
            key_schema = build_schema(val)
            properties[key] = {"type": "object", "properties": key_schema}
        else:
            print(f"Cannot determine type of {key} default input, {val}")
            pass # if type is not captured here, leave item out of the schema
    return properties

def build_rule_params_schema(full_schema, params_schema_template, paramsfile):
    params = read_yaml(paramsfile)
    program_names = list(params.keys())
    for name in program_names:
        prog_schema = params_schema_template['properties']['__prog__'].copy()
        prog_schema['properties']['program_params']['properties'] = build_schema(params[name]['program_params'])
        full_schema['properties'][name] = prog_schema
    return full_schema

def build_params_schema(args):
   # first, let's read in the pipeline defaults schema
    schema_dir = os.path.join(elvers_dir, "schemas")
    schema = read_yaml(find_input_file(args.defaults_template, name="pipeline defaults schema", add_paths=[schema_dir]))
    # the elvers params adds to the required sections of pipeline defaults. read in & update schema
    elvers_extra = read_yaml(find_input_file(args.elvers_template, name="elvers params schema", add_paths=[schema_dir]))
    update_nested_dict(schema, elvers_extra)
     # ok, now let's build schema for all
    rule_template = read_yaml(find_input_file(args.rule_template, name="rule params schema template", add_paths=[schema_dir]))
    # how do locations work with installed packages!?
    rule_paramsfiles = glob.glob(os.path.join(elvers_dir, 'rules', '*','params.yml'))
    for paramsfile in rule_paramsfiles:
        schema = build_rule_params_schema(schema, rule_template, paramsfile)

    # some tweaking using current targets. Set workflows to the current workflows. Set required to the rules in use
    targs = args.targets  #['get_data', 'trinity']

    # this is not exactly what's ending up in there. Not sure we need this.
    workflow_properties={'elvers':{'targets':targs}, 'required': ['elvers']}

    schema['properties']['elvers_workflows'] = {'type': 'object', 'properties': workflow_properties}#{'elvers: {'targets': targs}, 'required': ['elvers']}}
    schema['required'] = schema['required'] + targs

    write_yaml(schema, args.schema_output)

if __name__ == '__main__':
    """ very simple attempt at auto generating yaml schema for elvers """
    p = argparse.ArgumentParser()
    p.add_argument('paramsfile')
    p.add_argument('--rule_template', default= "rule_params.schema.yaml")
    p.add_argument('--elvers_template', default= "elvers_params.schema.yaml")
    p.add_argument('--defaults_template', default= "pipeline_defaults.schema.yaml")
    p.add_argument('--targets', nargs='*', default=[])
    p.add_argument('-o', '--schema_output', default = os.path.join(elvers_dir, 'schemas', 'elvers.fullschema.yaml'))
    args = p.parse_args()
    sys.exit(build_params_schema(args))
