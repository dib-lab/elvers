from os.path import join

def get_targets(units, assembly_basename, outdir, extensions = ['']):
    return [join(outdir, 'run_busco_' + assembly_basename)] #, "full_table_txome_busco.tsv")]
    #return [join('run_busco_' + assembly_basename, "full_table_txome_busco.tsv")]


