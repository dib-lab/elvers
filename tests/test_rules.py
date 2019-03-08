import subprocess
import os
import glob
import tempfile
import shutil
import pytest
import yaml

from utils import TempDirectory
from ..ep_utils import read_yaml, write_yaml,  update_nested_dict

def run_snakemake(homedir, tempdir, testdata, input_params, cmd, short=True, tests_dir = 'tests', rules_dir = 'rules'):
    # homedir is main eelpond dir, where conda envs should live
    conda_prefix = os.path.join(homedir, '.snakemake') 
    #testdata='salmon/test/assembly'
    shutil.copytree(os.path.join(homedir, rules_dir, testdata), os.path.join(tempdir, testdata))
    paramsfile = os.path.join(homedir, rules_dir, 'params.yml')
    default_params = read_yaml(paramsfile)

    ### now we need to read in and manipulate params





    snake_rule = os.path.realpath(os.path.join(rules_dir, rule))
    # specify fullpath to snakemake rule file
    cmd = cmd.extend(['-s', snake_rule])
    # short tests just do dryrun
    if short:
        cmd = cmd.append('-n')
    os.chdir(tempdir)
    subprocess.check_call(cmd) # run the command!
    # change back to main test dir
    os.chdir(origdir)
    

# this needs to change 

def run(ruledir, cmd, check_log=None):
    thisdir = os.getcwd() # we shoudl be in the tests dir.



    # main rules dir, all rules should be in here.
    rulesdir = os.path.join(os.dirname(thisdir), 'rules')
    # find everything required for tests
    rule = glob.glob(os.path.join(ruledir, '*.rule'))  # should only be one rule, but add a check for multiple?
    env = os.path.join(ruledir, 'environment.yaml')
    params = os.path.join(ruledir, 'params.yaml')
    testdir = os.path.join(ruledir, "test")

  # first, let's just try doing this within the testdir. After working, switch to a tempdir
    os.chdir(testdir)
    if os.path.exists(os.path.join(testdir,".snakemake")):
        shutil.rmtree(os.path.join(testdir,".snakemake"))
    subprocess.check_call(["snakemake", "--version"])
    try:
        subprocess.check_call(cmd) # run the command!
    except Exception as e:
       # go back to original directory
        os.chdir(origdir)
    
        # now grab logs
        logfiles = [os.path.join(d, f)
        for d, _, files in os.walk(os.path.join(testdir, "logs"))
            for f in files]
        for path in logfiles:
            with open(path) as f:
                msg = "###### Logfile: " + path + " ######"
                print(msg, "\n")
                print(f.read())
                print("#" * len(msg))
        if check_log is not None:
            for f in logfiles:
                check_log(open(f).read())
        else:
            raise e
    finally:
            # go back to original directory
        os.chdir(origdir)
    

# short (-n) vs long tests done

# better tempdir setup
#with utils.TempDirectory() as location:
#        testdata1 = utils.get_test_data('short.fa')
#        testdata2 = utils.get_test_data('short2.fa')
#        testdata3 = utils.get_test_data('short3.fa')
#        sigfile = os.path.join(location, 'short.fa.sig')
#
#        status, out, err = utils.runscript('sourmash',
#                                           ['compute', '-k', '31', '-o', sigfile,
#                                            testdata1,
#                                            testdata2, testdata3],
#                                           in_directory=location)


        #os.chdir(d) # I think we want to go into the ruledir, not testdir
       # run rule, testdata, 
        #dst = os.path.join(d, "master", wrapper)
        #    os.makedirs(dst, exist_ok=True)
        # need to set testdir as the outdir though! --> do defaults via test yaml?
       #cmd = cmd + ["--wrapper-prefix", "file://{}/".format(d)]
            
        # d was location of the tempdir --> we don't actually need it, unless we want to run tests elsewhere

#def test_trimmomatic_pe():
#    run("trimmomatic",
#        ["snakemake", "trimmed/a_1.fastq.gz", "--use-conda", "-F", "-s", "trimmomatic_pe.snake"])

#def test_trimmomatic_se():
#    run("trimmomatic",
#        ["snakemake", "trimmed/a.fastq.gz", "--use-conda", "-F", "-s", "trimmomatic_se.snake"])

def test_salmon_index():
    with TempDirectory() as location:
        # copy test data
        # add path to snakemake rule
        origdir = os.getcwd()
        conda_dir = os.path.join(origdir, '.snakemake')

        testdata='salmon/test/assembly'
        shutil.copytree(testdata, os.path.join(location, 'salmon/test/assembly'))
        snakefile = os.path.realpath(os.path.join(origdir, 'salmon', 'test', 'testing.snake'))
        os.chdir(location)
        run("salmon",
        ["snakemake", "salmon/transcriptome_index",  "--use-conda", "--conda-prefix", conda_dir, "-F", "-s", snakefile])
        os.chdir(origdir)


def test_salmon_quant_se():
    run("salmon",
        ["snakemake", "salmon/a_se_x_transcriptome/quant.sf",  "--use-conda", "-F", "-s", "testing.snake"])

def test_salmon_quant_pe():
    run("salmon",
        ["snakemake", "salmon/ab_pe_x_transcriptome/quant.sf",  "--use-conda", "-F", "-s", "testing.snake"])

