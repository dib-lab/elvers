import os
def get_targets(units, outdir):
    """
    Use the sample info provided in the tsv file
    to generate required targets for each workflow
    """
    #pre_trim_fastqc_targets = []
    post_trim_fastqc_targets = []
    for s, u in units.iterrows():
        sample, unit, read_type = u['sample'],u['unit'],u['read_type']
        post_trim_fastqc_targets.append(os.path.join(outdir, '{}_{}_1.trim_fastqc.zip'.format(sample,unit)))
        post_trim_fastqc_targets.append(os.path.join(outdir, '{}_{}_1.trim_fastqc.html'.format(sample,unit)))
        if read_type == 'pe':
            post_trim_fastqc_targets.append(os.path.join(outdir, '{}_{}_2.trim_fastqc.zip'.format(sample,unit)))
            post_trim_fastqc_targets.append(os.path.join(outdir, '{}_{}_2.trim_fastqc.html'.format(sample,unit)))
    return post_trim_fastqc_targets
