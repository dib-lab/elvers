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

def write_schema(yamlD, paramsfile):
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

def build_params_schema(paramsfile, outfile, rules = [], targets = [], rule_template =None, elvers_schema = None, pipeline_schema = None):
   # some setup, so we can use this within elvers main
    schema_dir = os.path.join(elvers_dir, "schemas")
    if not rule_template:
        rule_template = "rule_params.schema.yaml"
    if not elvers_schema:
        elvers_schema = "elvers_params.schema.yaml"
    if not pipeline_schema:
        pipeline_schema = "pipeline_defaults.schema.yaml"

   # first, let's read in the pipeline defaults schema
    schema = read_yaml(find_input_file(pipeline_schema, name="pipeline defaults schema", add_paths=[schema_dir]))
    # the elvers params adds to the required sections of pipeline defaults. read in & update schema
    elvers_extra = read_yaml(find_input_file(elvers_schema, name="elvers params schema", add_paths=[schema_dir]))
    update_nested_dict(schema, elvers_extra)
     # ok, now let's build schema for all
    for rule in rules:
        if rule in ['get_data', 'get_reference']:
            rule = 'utils'
        # getting carryover between rules. to avoid, read in fresh each time.
        r_template = read_yaml(find_input_file(rule_template, name="rule params schema template", add_paths=[schema_dir]))
        paramsfile = glob.glob(os.path.join(elvers_dir, 'rules', rule,'params.yml'))[0]
        schema = build_rule_params_schema(schema, r_template, paramsfile)

    # this is not exactly what we need, but it gets rid of `default` as a required target, so leave in for now.
    workflow_properties={'elvers':{'targets':targets}, 'required': ['elvers']}
    schema['properties']['elvers_workflows'] = {'type': 'object', 'properties': workflow_properties}
    schema['required'] = schema['required'] + rules

    write_schema(schema, outfile)

if __name__ == '__main__':
    """ very simple attempt at auto generating yaml schema for elvers """
    p = argparse.ArgumentParser()
    p.add_argument('paramsfile')
    p.add_argument('--rule_template', default= "rule_params.schema.yaml")
    p.add_argument('--elvers_template', default= "elvers_params.schema.yaml")
    p.add_argument('--defaults_template', default= "pipeline_defaults.schema.yaml")
    p.add_argument('--rules', nargs='*', default=[])
    p.add_argument('--targets', nargs='*', default=[])
    p.add_argument('-o', '--schema_output', default = os.path.join(elvers_dir, 'schemas', 'elvers.fullschema.yaml'))
    args = p.parse_args()
    sys.exit(build_params_schema(args.paramsfile, args.schema_output, args.rules, args.targets, args.rule_template, args.elvers_template, args.defaults_template))
