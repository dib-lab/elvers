from os import path

def get_targets(assembly_basename, outdir, extensions = ['.sig']):
    """
    Use the assembly_basename 
    to generate required targets for sourmash
    """
    sourmash_targs = []
    sourmash_targs = [assembly_basename + i for i in extensions] 
    return [path.join(outdir, targ) for targ in sourmash_targs]
