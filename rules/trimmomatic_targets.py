from os import path 

def get_targets(units, basename, outdir, extensions = ['.trim.fq.gz'], se_ext = ['se'], pe_ext = ['1','2']):
    """
    Use the sample info provided in the tsv file
    to generate required targets for trimmomatic
    """
    trim_targs = []
    for s, u in units.iterrows():
        sample, unit, read_type = u['sample'],u['unit'],u['read_type']
        end = se_ext if read_type == 'se' else pe_ext
        trim_targs = trim_targs +  ['{}_{}_'.format(sample,unit) + i + j for i in end for j in extensions]

    return [path.join(outdir, targ) for targ in trim_targs]

