from os import path
from common.utils import is_se

def get_targets(units, assembly_basename, outdir, extensions = ['.bam'], se_ext = ['se'], pe_ext = ['pe']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for bowtie2
    """
    bt2_targs = []

    # here we need to get all units belonging to a single sample.
    # When I used 'groupby' on our initial test data, groups were
    # A/AB, B/AB, rather than A, B, AB. I'm sure there's a way to
    # do exact matching in groupby, but here's a hackaround for now.
    samples = list(set(units['sample'].tolist()))
    #by_sample = units.groupby(['sample'], sort = False)
    for s in samples:
        unit_list = units.groupby(level=0).get_group(s)['unit'].tolist()
        for unit in unit_list:
            if is_se(units,s,unit):
                bt2_targs = bt2_targs + ['{}_{}_x_{}'.format(s, se_ext[0], assembly_basename) + i for i in extensions]
            else:
                bt2_targs = bt2_targs + ['{}_{}_x_{}'.format(s, pe_ext[0], assembly_basename) + i for i in extensions]
    bt2_targs = list(set(bt2_targs)) # elim any redundant targs
    bt2_index_targs = [assembly_basename + '_bowtie2' + i for i in ['.1.bt2','.rev.1.bt2']]
    bt2_tags = bt2_targs + bt2_index_targs
    return [path.join(outdir, targ) for targ in bt2_targs]
