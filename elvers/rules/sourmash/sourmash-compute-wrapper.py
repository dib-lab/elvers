"""Snakemake wrapper for sourmash compute."""

__author__ = "Lisa K. Johnson and N. Tessa Pierce"
__copyright__ = "Copyright 2019, Lisa K. Johnson and N. Tessa Pierce"
__email__ = "ljcohen@ucdavis.edu and ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
scaled = snakemake.params.get("scaled", "2000")
k = snakemake.params.get("k", "31")

k = [k] if (isinstance(k, str) or isinstance(k, int)) else k
k = ",".join(map(str, k))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "sourmash compute --scaled {scaled} -k {k} {snakemake.input} -o {snakemake.output} -p {snakemake.threads} {extra} {log}"
)
