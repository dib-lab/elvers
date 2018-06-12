from os import path

def get_targets(units, assembly_basename, outdir, extensions = ['/quant.sf', '/lib_format_counts.json']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    salmon_targs = []
    for s, u in units.iterrows():
        sample, unit, read_type = u['sample'],u['unit'],u['read_type']
        salmon_targs = salmon_targs +  ['{}_{}_x_{}'.format(sample,unit, assembly_basename) + i for i in extensions]

    return [path.join(outdir, targ) for targ in salmon_targs]


#or:
def get_left(wildcards):
    all_trimmed = expand(join(TRIM_DIR, '{sample}_{unit}_{end}.fq.gz'),**wildcards)
    right = [x for x in all_trimmed if '_2.trim' in x]
    left = [x.replace('_2.trim', '_1.trim') for x in right]
    single = [x for x in all_trimmed if '_1.trim' in x if x not in left]
    return right

def get_right(wildcards):
    all_trimmed = expand(join(TRIM_DIR, '{sample}_{unit}_{end}.fq.gz'),**wildcards)
    right = [x for x in all_trimmed if '_2.trim' in x]
    left = [x.replace('_2.trim', '_1.trim') for x in right]
    single = [x for x in all_trimmed if '_1.trim' in x if x not in left]
    return right

def get_single(wildcards):
    all_trimmed = expand(join(TRIM_DIR, '{sample}_{unit}_{end}.fq.gz'),**wildcards)
    right = [x for x in all_trimmed if '_2.trim' in x]
    left = [x.replace('_2.trim', '_1.trim') for x in right]
    single = [x for x in all_trimmed if '_1.trim' in x if x not in left]
    return single


