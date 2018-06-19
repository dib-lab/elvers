from os import path
from common.utils import is_se

def get_targets(units, basename, outdir, extensions = ['_fastqc.zip','_fastqc.html','.trim_fastqc.zip','.trim_fastqc.html'], se_ext = ['se'], pe_ext = ['1','2']):

    """
    Use the sample info provided in the tsv file
    to generate required targets for each workflow
    """
    fastqc_targs = []
    for s, u in units.iterrows():
        sample, unit = u['sample'],u['unit']
        end = se_ext if is_se(units,sample,unit) else pe_ext
        fastqc_targs = fastqc_targs +  ['{}_{}_'.format(sample,unit) + i + j for i in end for j in extensions]
    return [path.join(outdir, targ) for targ in fastqc_targs]
