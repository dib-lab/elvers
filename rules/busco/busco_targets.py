from os import path

def get_targets(units, assembly_basename, outdir, extensions = ['']):
    return [join(outdir, 'busco', 'run_busco_' + assembly_basename, "full_table_txome_busco.tsv")]


