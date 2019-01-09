__author__ = "N .Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
outdir = path.dirname(snakemake.output.get('r1'))

r1 = snakemake.input.get("r1")
r2 =  snakemake.input.get("r2")
r = snakemake.input.get("r")

def move_files(outdir, in_list, out_list):
    for f, o in zip(in_list, out_list):
        print(f)
        print(o)
        f = path.join(outdir, f)
        o = path.join(outdir, o)
        shell("cp {f} {o}")
        shell("rm -f {f}")
    
def build_default_outname(infile):
    # Rcorrector outputs gzipped files IF input files are gzipped
    end = '.gz' if infile.endswith('.gz') else ''
    return(path.basename(infile.rsplit('.f')[0]) + '.cor.fq' + end)

assert (r1 is not None and r2 is not None) or r is not None, "either r1 and r2 (paired), or r (unpaired) are required as input"
if r1:
    # handle inputs
    r1 = [snakemake.input.r1] if isinstance(snakemake.input.r1, str) else snakemake.input.r1
    r2 = [snakemake.input.r2] if isinstance(snakemake.input.r2, str) else snakemake.input.r2
    assert len(r1) == len(r2), "input-> equal number of files required for r1 and r2"
    r1_cmd = ' -1 ' + ",".join(r1)
    r2_cmd = ' -2 ' + ",".join(r2)
    read_cmd = " ".join([r1_cmd,r2_cmd])
    # handle outputs
    r1_out = [snakemake.output.r1] if isinstance(snakemake.output.r1, str) else snakemake.output.r1
    r2_out = [snakemake.output.r2] if isinstance(snakemake.output.r2, str) else snakemake.output.r2
    r1_default, r2_default = [], []
    for f in r1:
        r1_default+= [build_default_outname(f)]
    for f in r2:
        r2_default+= [build_default_outname(f)]
if r:
    # handle inputs
    assert r1 is None and r2 is None, "cannot handle mixed paired/unpaired input files. Please input either r1,r2 (paired) or r (unpaired)"
    r = [snakemake.input.r] if isinstance(snakemake.input.r, str) else snakemake.input.r
    read_cmd = ' -r ' + ",".join(r)
    # handle outputs
    r_out = [snakemake.output.r] if isinstance(snakemake.output.r, str) else snakemake.output.r
    r_default = []
    for f in r:
        r_default += [build_default_outname(f)]

shell("run_rcorrector.pl {read_cmd} -od {outdir} {snakemake.params.extra} -t {snakemake.threads} {log}")

if r1_default:
    move_files(outdir, r1_default, r1_out)
    move_files(outdir, r2_default, r2_out)
elif r_default:
    move_files(outdir, r_default, r_out)
