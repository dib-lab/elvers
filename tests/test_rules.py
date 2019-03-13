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
from .const import (here, test_config_yaml)

def run_ruletest(rulename, testdir, extra_configD = {}, short = True): # can we pass in rulename, paramsD here, testdata, short yes/no? 
    """ test a rule or workflow"""
    # set up dirs
    homedir = os.path.dirname(here)
    run_eelpond_cmd = os.path.join(homedir, 'run_eelpond')
    conda_prefix = os.path.join(homedir, '.snakemake')
    
    # test info from rule
    rulefile = glob.glob(os.path.join(homedir, 'rules', '*', rulename + '.rule'))[0]
    ruledir = os.path.dirname(rulefile)
    test_yml = glob.glob(os.path.join(ruledir, testdir,'*.yml'))[0]

    with TempDirectory() as location:
        # copy in test data
        #tmpdir = os.path.join(location, testdir)
        #shutil.copytree(os.path.join(ruledir, testdir), tmpdir)
        #test_yml = glob.glob(os.path.join(tmpdir, '*.yml'))[0]
        
        ## or, avoid copying data this way: ##
        os.chdir(os.path.join(ruledir, testdir))
        # need to be here in order to properly find any relative assemblyinput, etc paths. Maybe fix this in run_eelpond to get path relative to file
        cmd = [run_eelpond_cmd, test_yml, rulename, 'get_data', '--conda_prefix', conda_prefix, '--out_path', location]
        ## NOTE: get data should be added automatically if we need the reads, via high-level input checks. Take out of here once that is implemented issue #110
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
            # alternatively use this?
            #p_out, p_err = capture_stdouterr(command,here)
        except Exception as e:
            # need to check and print the logs here?
            raise e
        finally:
            os.chdir(here) # back to tests dir


def test_salmon():
     run_ruletest('salmon', 'test', {})
     run_ruletest('salmon', "test", {'salmon':{'program_params': {'quant_params':{'libtype': "IU"}}}})

def test_salmon_long():
     run_ruletest('salmon', 'test', {}, short=False)



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

