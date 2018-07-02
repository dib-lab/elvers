from os import path
from common.utils import is_se

def get_targets(samples, basename, outdir, extensions = ['.diffexp.tsv', '.ma-plot.svg', ], se_ext = ['se'], pe_ext = ['1','2']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    de_targs = []
    for sample, contrast in samples.iterrows():
        sample, contrast = u['sample'],u['contrast']
        de_targs = de_targs +  ['{}'.format(contrast) + i for i in extensions]
    de_targs = de_targs + ['pca.svg']
    return [path.join(outdir, targ) for targ in de_targs]

# targs from rna-seq-star
# expand(["results/diffexp/{contrast}.diffexp.tsv",
#                "results/diffexp/{contrast}.ma-plot.svg"],
#               contrast=config["diffexp"]["contrasts"]),
#        "results/pca.svg"



