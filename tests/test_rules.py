 #!/usr/bin/env python
import subprocess
import os
import glob
import tempfile
import shutil
import pytest
import yaml

from .utils import TempDirectory
from .const import (here, run_eelpond_cmd, test_config_yaml)
from .ep_utils import read_yaml, write_yaml, find_yaml, update_nested_dict

def grab_pipeline_defaults(homedir):
    pipeline_defaultsFile = find_yaml(homedir, os.path.join('ep_utils', 'pipeline_defaults'), 'pipeline_defaults')
    pipeline_defaults = read_yaml(pipeline_defaultsFile)
    return pipeline_defaults

def run_ruletest(rulename, testdata, rule_output, input_configD = {}, short = True): # can we pass in rulename, paramsD here, testdata, short yes/no? 
    """Test a snakemake rule"""
    # homedir is main eelpond dir, where conda envs should live
    homedir = os.path.dirname(here)
    conda_prefix = os.path.join(homedir, '.snakemake')
    rulesdir = os.path.join(homedir, 'rules')
    cmd = ["snakemake", rule_output, "--conda_prefix", conda_prefix, "--use-conda"]
    # short tests just do dryrun
    if short:
        cmd.append('-n')
    # ^ above here is the same for every rule test --> make setup for test class?
    
    # SETUP FOR THE SPECIFIC RULE --> into separate function?
    rulefile = glob.glob(os.path.join(rulesdir, '*', rulename + '.rule'))[0]
    cmd = cmd.extend(['-s', rulefile]) # specify fullpath to snakemake rule file 
    ruledir = os.path.dirname(rulefile)
    # grab default directory structure, etc
    paramsD = grab_pipeline_defaults(homedir) 
    program_defaultsfile = os.path.join(ruledir, 'params.yml')
    import pdb;pdb.set_trace()
    paramsD[rulename] = read_yaml(program_defaultsfile).get(rulename) # opens us up to multiple programs per params.yml, as likely for utils
    
    # now update default params using input_paramsD
    if input_configD: # need to make sure this is appropriately formatted. Maybe do some manipulation here.
        update_nested_dict(paramsD, input_configD)

    with TempDirectory() as location:
        # copy testdata to tmpdir
        shutil.copytree(os.path.join(rules_dir, testdata), os.path.join(location, testdata))
        # print a new params file and add this to the snakemake command
        test_paramsfile = os.path.join(location,'config.yml')
        
        write_yaml(test_paramsfile, paramsD)
        cmd = cmd.extend(['--configfile', test_paramsfile]) # specify fullpath to snakemake rule file 
        
        print(cmd)
        
        # Actually run the test
        # change into temp dir and run!
        os.chdir(location)
        subprocess.check_call(["snakemake", "--version"]) # inside the try block?
        try:
            subprocess.check_call(cmd) # run the command!
            # alternatively use this?
            #p_out, p_err = capture_stdouterr(command,here)
        except Exception as e:
            # want all logs to be right in here. set it up earlier.
            #logfiles = [os.path.join(d, f) for d, _, files in os.walk("logs") for f in files]
            #for path in logfiles:
            #    with open(path) as f:
            #        msg = "###### Logfile: " + path + " ######"
            #        print(msg, "\n")
            #        print(f.read())
            #        print("#" * len(msg))
            #for f in logfiles:
            #    check_log(open(f).read())
            #else:
            raise e
        finally:
            os.chdir(here) # back to tests dir


input_params = {}

def test_salmon_index():
     run_ruletest('salmon', "salmon/transcriptome_index", 'salmon/test/assembly', {})

#def test_salmon_index():
#    with TempDirectory() as location:
        # copy test data
        # add path to snakemake rule
#        origdir = os.getcwd()
#        conda_dir = os.path.join(origdir, '.snakemake')

#        testdata='salmon/test/assembly'
#        shutil.copytree(testdata, os.path.join(location, 'salmon/test/assembly'))
#        snakefile = os.path.realpath(os.path.join(origdir, 'salmon', 'test', 'testing.snake'))
#        os.chdir(location)
#        run("salmon",
#        ["snakemake", "salmon/transcriptome_index",  "--use-conda", "--conda-prefix", conda_dir, "-F", "-s", snakefile])
#        os.chdir(origdir)


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

