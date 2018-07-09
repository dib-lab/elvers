from os.path import join

def get_targets(units, basename, outdir, extensions = ['.fasta', '.fasta.gene_trans_map']):
    """
    Use run basename from config
    to generate Trinity targets
    """
    mv_assemb_targs = ['{}'.format(basename) + i for i in extensions]
    return [join(outdir, targ) for targ in mv_assemb_targs]


