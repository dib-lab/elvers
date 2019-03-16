 #!/usr/bin/env python
import subprocess
import os
import glob
import tempfile
import shutil
import pytest
import yaml
import pandas as pd

from .utils import TempDirectory
from .const import (here, run_eelpond_cmd, test_config_yaml)
from .ep_utils import * 

def grab_pipeline_defaults(homedir):
    pipeline_defaultsFile = find_yaml(homedir, os.path.join('ep_utils', 'pipeline_defaults'), 'pipeline_defaults')
    pipeline_defaults = read_yaml(pipeline_defaultsFile)
    return pipeline_defaults

def run_ruletest(rulename, test_yml, testdata, extra_configD = {}, short = True): # can we pass in rulename, paramsD here, testdata, short yes/no? 
    """Test a snakemake rule. Basically a lightweight version of run_eelpond, for testing single rules"""
    # homedir is main eelpond dir, where conda envs should live
    homedir = os.path.dirname(here)
    conda_prefix = os.path.join(homedir, '.snakemake')
    rulesdir = os.path.join(homedir, 'rules')
    # ^ above here is the same for every rule test --> make setup for test class?
    
    # SETUP FOR THE SPECIFIC RULE --> into separate function?
    rulefile = glob.glob(os.path.join(rulesdir, '*', rulename + '.rule'))[0]
    ruledir = os.path.dirname(rulefile)
    testdir = os.path.join(ruledir, 'test')
    test_config = read_yaml(os.path.join(testdir, test_yml))
    
    # grab default directory structure, etc
    paramsD = grab_pipeline_defaults(homedir) 
    paramsD['eelpond_workflows'] = {rulename:{'targets':[rulename]}} # if we want to test workflows too, change this.
    program_defaultsfile = os.path.join(ruledir, 'params.yml')
    paramsD[rulename] = read_yaml(program_defaultsfile).get(rulename) # opens us up to multiple programs per params.yml, as likely for utils

    # read in test yml (specifies samples, assemblyinput if needed)
    os.chdir(testdir)
    paramsD.update(test_config)
    samples = pd.read_csv(paramsD["samples"],dtype=str, sep='\t').set_index(["sample", "unit"], drop=False)
    samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
    #validate(samples, schema="schemas/samples_v2.schema.yaml") # new version
    assemblyinput = paramsD.get('assemblyinput', None)
    if assemblyinput:
        paramsD['assemblyinput'] = {'program_params': assemblyinput} 
        paramsD, assemb_ext = handle_assemblyinput(assemblyinput, paramsD)
        paramsD['assembly_extensions'] = list(set(paramsD.get('assembly_extensions', []) + [assemb_ext]))
        paramsD['eelpond_workflows'][rulename]['targets'].append('assemblyinput') 

    # now update default params using input_paramsD
    if extra_configD: # need to make sure this is appropriately formatted. Maybe do some manipulation here.
        update_nested_dict(paramsD, extra_configD)
    rule_outputs = generate_all_targs(paramsD, samples)
    print(rule_outputs)
    rule_output = ','.join(rule_outputs)
    cmd = ["snakemake", rule_outputs[0], "--conda-prefix", conda_prefix, "--use-conda"]
    cmd.extend(['-s', rulefile]) # specify fullpath to snakemake rule file 
    # short tests just do dryrun
    if short:
        cmd.append('-n')
    
    with TempDirectory() as location:
        # copy testdata to tmpdir
        shutil.copyfile(os.path.join(testdir, paramsD['samples']),os.path.join(location, paramsD['samples']))
        shutil.copytree(os.path.join(testdir, 'assembly'), os.path.join(location, 'assembly'))
        shutil.copytree(os.path.join(testdir, 'reads'), os.path.join(location, 'reads'))
        os.chdir(location)
        # print a new params file and add this to the snakemake command
        test_paramsfile = os.path.join(location,'config.yml')
        write_yaml(paramsD, test_paramsfile)
        import pdb;pdb.set_trace()
        cmd.extend(['--configfile', test_paramsfile]) # specify fullpath to snakemake rule file 
        
        # change into temp dir and actually run!
        os.chdir(location)
        subprocess.check_call(["snakemake", "--version"]) # inside the try block?
        try:
            subprocess.check_call(cmd) # run the command!
            # alternatively use this?
            #p_out, p_err = capture_stdouterr(command,here)
        except Exception as e:
            # need to check and print the logs here?
            raise e
        finally:
            os.chdir(here) # back to tests dir


def test_salmon_index():
     run_ruletest('salmon', 'test.yml', 'assembly', {})
     #run_ruletest('salmon', "quant/transcriptome.salmonindex", 'assembly', {'salmon':{'program_params': {'quant_params':{'libtype': "IU"}}}})

#def test_salmon_quant_se():
#    run("salmon",
#        ["snakemake", "salmon/a_se_x_transcriptome/quant.sf",  "--use-conda", "-F", "-s", "testing.snake"])

#def test_salmon_quant_pe():
#    run("salmon",
#        ["snakemake", "salmon/ab_pe_x_transcriptome/quant.sf",  "--use-conda", "-F", "-s", "testing.snake"])

# test creating and updating environments
#def run_test_env(rulename):
#    cmd = cmd.append('--create_envs_only')
#    
#    if os.path.exists(os.path.join(testdir,".snakemake")):
#        shutil.rmtree(os.path.join(testdir,".snakemake"))
#    subprocess.check_call(["snakemake", "--version"])
#    try:
#        subprocess.check_call(cmd) # run the command!
#    except Exception as e:
       # go back to original directory

