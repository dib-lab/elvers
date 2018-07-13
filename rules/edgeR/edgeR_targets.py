from os.path import join
from common.utils import get_params


def get_targets(units, basename, outdir, conf = '', extensions = ['.diffexp.tsv', '.ma-plot.pdf', ], se_ext = ['se'], pe_ext = ['1','2']):
    """
    generate contrasts from deseq2 configfile #not best soln - think more.
    """
    # default contrast info from default params
    contrast_info = get_params('edgeR')
    # override defaults with snakemake config (main configfile or snakemake --config)
    # not gonna work - haven't gotten config yet in main snakefile?
    if 'edgeR' in conf.keys():
        contrast_info.update(conf['edgeR'])
#    import pdb;pdb.set_trace()
    contrast_list = contrast_info['contrasts']
    # build targets
    de_targs = []
    for c in contrast_list:
        de_targs = de_targs +  ['{}'.format(c) + i for i in extensions]
    de_targs = de_targs + ['all.rds']
    de_targs = de_targs + ['pca.pdf']
    #de_targs = de_targs + ['pca.svg']
    return [join(outdir, targ) for targ in de_targs]
