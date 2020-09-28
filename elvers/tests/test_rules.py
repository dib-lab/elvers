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

elvers_dir = os.path.dirname(os.path.dirname(here))


def run_ruletest(rulename, extra_configD = {}, protein_ref = False, short = True, envs_only=False): # can we pass in rulename, paramsD here, testdata, short yes/no?
    """ test a rule or workflow"""
    # set up dirs
    conda_prefix = os.path.join(elvers_dir, '.snakemake')
    # test info from rule
    rulefile = glob.glob(os.path.join(elvers_dir, 'elvers/rules', '*', rulename + '.rule'))[0]
    ruledir = os.path.dirname(rulefile)
    testdir = os.path.join(here, "test_files")
    if protein_ref:
        test_yml = os.path.join(testdir, 'prot_test.yml')
    elif short:
        test_yml = os.path.join(testdir, 'short_test.yml')
    else:
        test_yml = os.path.join(testdir, 'long_test.yml')
    try:
        additional_test_yml = glob.glob(os.path.join(ruledir, 'test','test.yml'))[0]
        add_params = ['--extra_config', additional_test_yml]
    except:
        additional_test_yml = None
        add_params = []

    # test creating and updating environments
    if envs_only:
        cmd = cmd.append('--create_envs_only')
   #    if os.path.exists(os.path.join(testdir,".snakemake")):
   #        shutil.rmtree(os.path.join(testdir,".snakemake"))

    with TempDirectory() as location:
        # copy in test data
        os.chdir(testdir)
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

def test_get_data_short():
    run_ruletest('get_data')

@pytest.mark.long
def test_get_data_long():
    run_ruletest('get_data', short = False)

def test_get_reference_short():
    run_ruletest('get_reference')

@pytest.mark.long
def test_get_reference_long():
    run_ruletest('get_reference', short = False)

def test_dammit_short():
    run_ruletest('dammit')

@pytest.mark.long
def test_dammit_long():
    db_dir = os.path.join(elvers_dir, 'databases')
    run_ruletest('dammit', extra_configD = {'dammit':{'db_dir': db_dir}}, short=False)

def test_salmon_short():
    run_ruletest('salmon')
    run_ruletest('salmon', extra_configD = {'salmon':{'quant_params':{'libtype': "IU"}}})

@pytest.mark.long
def test_salmon_long():
     run_ruletest('salmon', short=False)

def test_trimmomatic_short():
    run_ruletest('trimmomatic')

@pytest.mark.long
def test_trimmomatic_long():
    run_ruletest('trimmomatic', short=False)

def test_khmer_short():
    run_ruletest('khmer')

@pytest.mark.long
def test_khmer_long():
    run_ruletest('khmer', short=False)

def test_trinity_short():
    run_ruletest('trinity')

@pytest.mark.long
def test_trinity_long():
    run_ruletest('trinity', short=False)

def test_plass_short():
    run_ruletest('plass')

@pytest.mark.long
def test_plass_long():
    run_ruletest('plass', short=False)

def test_pear_short():
    run_ruletest('pear')

@pytest.mark.long
def test_pear_long():
    run_ruletest('pear', short = False)

def test_paladin_short():
    run_ruletest('paladin', protein_ref = True)

@pytest.mark.long
def test_paladin_long():
    run_ruletest('paladin', protein_ref = True, short = False)

def test_rcorrector_short():
    run_ruletest('rcorrector')

@pytest.mark.long
def test_rcorrector_long():
    run_ruletest('rcorrector', short = False)

def test_bowtie2_short():
    run_ruletest('bowtie2')

@pytest.mark.long
def test_bowtie2_long():
    run_ruletest('bowtie2', short = False)


