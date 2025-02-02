import os
import glob
import pandas as pd
import argparse
from collections import defaultdict


"""
This should be ran AFTER alignment to the reference sequences where the nt length of the reference is maintained.

This code will: 
1. first check length of reference
2. Check length of all other sequences in dictionary
3. Remove the 5' and 3' UTR from all sequences (including ref)

python3 adjust_haplotype_names.py --fastas ../../projects/autopsy_chertow/haplotypes/fasta10/ --outfile cat.haplotypes.10.fasta
"""
parser = argparse.ArgumentParser()
parser.add_argument('--fastas','-f',required=True,help='Indicate path and reference fasta file') #args.ref
parser.add_argument('--outfile',default = 'SARS.fasta', help='Indicate directory to save coverage csv') #args.ref
args = parser.parse_args()

def read_fasta(fp): #function to read each segment of fasta, from tim readreport
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def open_fasta(filename): #opens each fasta- makes a dictionary of segment name from tim readreport
    # filename = '../FILES/reference/'+someref
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict

outnt = open(args.outfile, 'w')
for filename in os.listdir(args.fastas):
    n = filename.split(".")[0]
    f = open_fasta("{0}/{1}".format(args.fastas, filename))
    for key, value in f.items():
        name = "{0}_{1}".format(n, key)
        outnt.write('>' + name + '\n' + value + '\n')

outnt.close()
