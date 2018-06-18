
from os import path 

def get_targets(units,basename, outdir, extensions = ['.fasta']):
    """
    Use run basename from config
    to generate Trinity targets
    """
    trinity_targs = ['{}'.format(basename) + i for i in extensions]

    return [path.join(outdir, targ) for targ in trinity_targs]

