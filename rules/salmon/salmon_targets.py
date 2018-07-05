from os import path

def get_targets(units, assembly_basename, outdir, extensions = ['/quant.sf', '/lib_format_counts.json'], se_ext = ['se'], pe_ext = ['pe']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for salmon
    """
    salmon_targs = []
    pe_only = units[units['fq2'].notnull()]
    for s, u in pe_only.iterrows():
        sample= u['sample']
        salmon_targs = salmon_targs + ['{}_{}_x_{}'.format(sample, pe_ext[0], assembly_basename) + i for i in extensions]
    se_only = units[units['fq2'].isnull()]
    for s, u in se_only.iterrows():
        sample= u['sample']
        salmon_targs = salmon_targs + ['{}_{}_x_{}'.format(sample, se_ext[0], assembly_basename) + i for i in extensions]
    
    salmon_targs = list(set(salmon_targs)) # elim any redundant targs
    return [path.join(outdir, targ) for targ in salmon_targs]
