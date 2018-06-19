from os.path import join
from common.utils import is_se

def get_targets(units, basename, outdir, extensions = ['.fasta']):
    """
    Use run basename from config
    to generate Trinity targets
    """
    trinity_targs = ['{}'.format(basename) + i for i in extensions]
    return [join(outdir, targ) for targ in trinity_targs]


def get_trimmed_trinity_input(units, basename, outdir, extensions = ['.trim.fq.gz'], se_ext = ['se'], pe_ext = ['1','2']):

    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    trim_targs = []
    for s, u in units.iterrows():
        sample, unit = u['sample'],u['unit']
        end = se_ext if is_se(units,sample, unit) else pe_ext
        trim_targs = trim_targs +  ['{}_{}_'.format(sample, unit) + i + j for i in end for j in extensions]
    return [join(outdir, targ) for targ in trim_targs]

