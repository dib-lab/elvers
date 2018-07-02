from os.path import join
from common.utils import get_params


def get_targets(samples, basename, outdir, extensions = ['.diffexp.tsv', '.ma-plot.svg', ], se_ext = ['se'], pe_ext = ['1','2']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    # default contrast info from deseq2 yaml 
    contrast_info = get_params('deseq2')
    # override defaults with snakemake config (main configfile or snakemake --config)
    if 'deseq2' in config.keys():
        contrast_info.update(config['deseq2'])
    contrast_list = contrast['contrasts']
    # build targets
    de_targs = []
    for c in contrast_list: 
        de_targs = de_targs +  ['{}'.format(c) + i for i in extensions]
    de_targs = de_targs + ['pca.svg']
    return [path.join(outdir, targ) for targ in de_targs]

# targs from rna-seq-star
# expand(["results/diffexp/{contrast}.diffexp.tsv",
#                "results/diffexp/{contrast}.ma-plot.svg"],
#               contrast=config["diffexp"]["contrasts"]),
#        "results/pca.svg"



