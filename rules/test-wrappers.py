import subprocess
import os
import glob
import tempfile
import shutil
import pytest
import yaml

# starting from https://bitbucket.org/snakemake/snakemake-wrappers/src/5a5bd45590896a7c7ca00bd6d558cdf40bc78c20/test.py?at=master&fileviewer=file-view-default

def run(ruledir, cmd, check_log=None):
    origdir = os.getcwd() # rules should all be subdirs in here
    
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
    
#    with tempfile.TemporaryDirectory() as d: 
        # copy test data and snakefile into tempdir
        #shutil.copy(env, d)
        #shutil.copy(params, d)
        
#        shutil.copytree(test, d)

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

class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='sourmashtest_')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


# to do: copy test data, change into it, run snakemake from there (specify path to snakefile)


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
    run("salmon",
        ["snakemake", "salmon/transcriptome_index",  "--use-conda", "-F", "-s", "testing.snake"])

def test_salmon_quant_se():
    run("salmon",
        ["snakemake", "salmon/a_se_x_transcriptome/quant.sf",  "--use-conda", "-F", "-s", "testing.snake"])

def test_salmon_quant_pe():
    run("salmon",
        ["snakemake", "salmon/ab_pe_x_transcriptome/quant.sf",  "--use-conda", "-F", "-s", "testing.snake"])

