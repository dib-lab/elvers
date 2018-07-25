from os.path import join
from common.utils import get_params


def get_targets(units, basename, outdir, conf = '', extensions = ['.diffexp.tsv', '.padj_0.1.diffexp.tsv', '.ma-plot.pdf', ], se_ext = ['se'], pe_ext = ['1','2']):
    """
    generate contrasts from deseq2 configfile #not best soln - think more.
    """
    # default contrast info from deseq2 yaml 
    contrast_info = get_params('deseq2')
    # override defaults with snakemake config (main configfile or snakemake --config)
    # not gonna work - haven't gotten config yet in main snakefile?
    if 'deseq2' in conf.keys():
        contrast_info.update(conf['deseq2'])
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


# targs from rna-seq-star
# expand(["results/diffexp/{contrast}.diffexp.tsv",
#                "results/diffexp/{contrast}.ma-plot.svg"],
#               contrast=config["diffexp"]["contrasts"]),
#        "results/pca.svg"



