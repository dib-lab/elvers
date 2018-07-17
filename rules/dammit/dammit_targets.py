from os import path
from common.utils import is_se

def get_targets(units, assembly_basename, outdir, extensions = ['.fasta.dammit.gff3'], se_ext = ['se'], pe_ext = ['pe']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for dammit
    """
    dammit_targs = [assembly_basename + i for i in extensions]
    return [path.join(outdir, targ) for targ in dammit_targs]
