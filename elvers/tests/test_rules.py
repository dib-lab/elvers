 #!/usr/bin/env python
import subprocess
import os
import glob
import tempfile
import shutil
import pytest
import yaml
import pandas as pd

from .utils import TempDirectory, write_yaml
from .const import (here, test_config_yaml, elvers_cmd)

def run_ruletest(rulename, testdir, extra_configD = {}, short = True): # can we pass in rulename, paramsD here, testdata, short yes/no?
    """ test a rule or workflow"""
    # set up dirs
    homedir = os.path.dirname(os.path.dirname(here))
    conda_prefix = os.path.join(homedir, '.snakemake')

    # test info from rule
    rulefile = glob.glob(os.path.join(homedir, 'elvers/rules', '*', rulename + '.rule'))[0]
    ruledir = os.path.dirname(rulefile)
    if short:
         test_yml = os.path.join(here, 'test_files', 'short_test.yml')
    else:
         test_yml = os.path.join(here, 'test_files', 'long_test.yml')

    try:
        additional_test_yml = glob.glob(os.path.join(ruledir, 'test','*.yml'))[0]
    except:
        additional_test_yml = None
    add_params = []
    if additional_test_yml:
        add_params = ['--extra_config', additional_test_yml]

    with TempDirectory() as location:
        # copy in test data
        os.chdir(os.path.join(here, testdir))
        cmd = [elvers_cmd, test_yml, rulename, '--conda_prefix', conda_prefix, '--out_path', location] + add_params
        # short tests just do dryrun
        if short:
            cmd.append('-n')
        if extra_configD:
            extra_yml = os.path.join(location, 'extra.yml')
            write_yaml(extra_configD, extra_yml)
            cmd.extend(['--extra_config', extra_yml])
            #cmd.extend(['--config_dict', extra_configD]) # not working via pytest
        subprocess.check_call(["snakemake", "--version"])
        try:
            subprocess.check_call(cmd) # run the command!
        except Exception as e:
            # need to check and print the logs here?
            raise e
        finally:
            os.chdir(here) # back to tests dir


def test_salmon_short():
     run_ruletest('salmon', 'test_files', {})
     run_ruletest('salmon', "test_files", {'salmon':{'program_params': {'quant_params':{'libtype': "IU"}}}})

def test_get_data_short():
     run_ruletest('get_data', 'test_files', {})

def test_get_reference_short():
     run_ruletest('get_reference', 'test_files', {})

def test_dammit_short():
     run_ruletest('dammit', 'test_files', {})

#def test_salmon_long():
     #run_ruletest('salmon', 'test', {}, short=False)



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

