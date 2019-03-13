import os
import argparse
import screed

def generate_geneTransMap(fasta, geneTransMap, transGeneMap):
    if transGeneMap:
        outT = open(transGeneMap, "w")
    with open(geneTransMap, "w") as outG:
        with screed.open(fasta) as seqs:
            for read in seqs:
                trans = read.name.split(' ')[0]
                gene = trans.rsplit('_', 1)[0]
                outG.write(gene + '\t' + trans + '\n')
                if transGeneMap:
                    outT.write(trans + '\t' + gene + '\n')

    if transGeneMap:
        outT.close()


if __name__ == '__main__':
    """Function: Take in a trinity assembly; output a gene to trans map. Optionally,
also output trans \t gene, for use with tximport of salmon files for gene-level DE.
"""
    psr = argparse.ArgumentParser()
    psr.add_argument('--fasta')
    psr.add_argument('--geneTransMap')
    psr.add_argument('--transGeneMap', default=None)
    args = psr.parse_args()
    generate_geneTransMap(args.fasta,args.geneTransMap,args.transGeneMap)