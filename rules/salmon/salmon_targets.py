from os import path
from common.utils import is_se

def get_targets(units, assembly_basename, outdir, extensions = ['/quant.sf', '/lib_format_counts.json'], se_ext = ['se'], pe_ext = ['pe']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for salmon
    """
    salmon_targs = []

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
                salmon_targs = salmon_targs + ['{}_{}_x_{}'.format(s, se_ext[0], assembly_basename) + i for i in extensions]
            else:
                salmon_targs = salmon_targs + ['{}_{}_x_{}'.format(s, pe_ext[0], assembly_basename) + i for i in extensions]
    salmon_targs = list(set(salmon_targs)) # elim any redundant targs
    return [path.join(outdir, targ) for targ in salmon_targs]
