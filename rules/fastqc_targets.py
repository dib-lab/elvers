from os import path

def get_targets(units, basename, outdir, extensions = ['_fastqc.zip','_fastqc.html','.trim_fastqc.zip','.trim_fastqc.html'], se_ext = ['1'], pe_ext = ['1','2']):

    """
    Use the sample info provided in the tsv file
    to generate required targets for each workflow
    """
    fastqc_targs = []
    for s, u in units.iterrows():
        sample, unit, read_type = u['sample'],u['unit'],u['read_type']
        end = se_ext if read_type == 'se' else pe_ext
        fastqc_targs = fastqc_targs +  ['{}_{}_'.format(sample,unit) + i + j for i in end for j in extensions]
    return [path.join(outdir, targ) for targ in fastqc_targs]
