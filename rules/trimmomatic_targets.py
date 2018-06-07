import os

def get_targets(units, TRIM_DIR): 
# here we could enable changes via config output format, if desired... 
    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    trim_targets = []
    for s, u in units.iterrows():
        sample, unit, read_type = u['sample'],u['unit'],u['read_type']
        trim_targets.append(os.path.join(TRIM_DIR, '{}_{}_1.trim.fq.gz'.format(sample,unit)))
        if read_type == 'pe':
            trim_targets.append(os.path.join(TRIM_DIR, '{}_{}_2.trim.fq.gz'.format(sample,unit)))
    return trim_targets
