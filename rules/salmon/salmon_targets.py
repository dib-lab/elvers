from os import path

def get_targets(units, assembly_basename, outdir, extensions = ['/quant.sf', '/lib_format_counts.json']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for salmon
    """
    salmon_targs = []
    for s, u in units.iterrows():
        sample, unit = u['sample'],u['unit']
        salmon_targs = salmon_targs +  ['{}_{}_x_{}'.format(sample,unit, assembly_basename) + i for i in extensions]

    return [path.join(outdir, targ) for targ in salmon_targs]
