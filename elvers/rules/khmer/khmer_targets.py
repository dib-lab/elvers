from os import path
from common.utils import is_se

def get_targets(units, basename, outdir, extensions = ['.gz'], se_ext = ['.se'], pe_ext = ['.paired', '.single', '.paired.1','.paired.2']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    targs = []
    for s, u in units.iterrows():
        sample, unit = u['sample'],u['unit']
        end = se_ext if is_se(units,sample,unit) else pe_ext
        targs = targs +  ['{}_{}_'.format(sample,unit) + i + j for i in end for j in extensions]
    return [path.join(outdir, targ) for targ in targs]
