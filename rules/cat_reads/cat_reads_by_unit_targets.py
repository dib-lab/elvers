from os import path
from common.utils import is_se

def get_targets(units, basename, outdir, extensions = ['.fq.gz'], se_ext = ['se'], pe_ext = ['1','2']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for cat_reads_by_unit
    """
    cat_targs = []
    for s, u in units.iterrows():
        sample, unit = u['sample'],u['unit']
        end = se_ext if is_se(units,sample,unit) else pe_ext
        cat_targs = cat_targs +  ['{}_'.format(sample) + i + j for i in end for j in extensions]

    return [path.join(outdir, targ) for targ in cat_targs]

