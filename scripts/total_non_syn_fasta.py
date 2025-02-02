"""
Written by: Kate Johnson


module purge
module load python/intel/3.8.6
module load pysam/0.16.0.1


# input: csv file with reference information only for coding regions. ORF1a/ORF1b (not nsps) in the case of SARS-CoV-2
# 
python3 total_non_syn_fasta.py --features 

"""
import os
import glob
import pandas as pd
import numpy as np
import math
import scipy
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--ref', required=True , help='flu ref')
parser.add_argument('--features', required=True, help='file to be used for variant stuff')
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
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict

def getindex(codon):
    """
    Identifying where the nucleotide falls in codon
    INPUT: major codon from snplist file
    OUTPUT: Index location of nucleotide
    """
    return [i for i, c in enumerate(codon) if c.isupper()]

aminoacid = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
            'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
            'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
            'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
            'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
            'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
            'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
            'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
            'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
            'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
            'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
            'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
            'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
            'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
            'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
            'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}

def adjust_codon(codon):
    codon_list = [codon] * 3 # repeat the codon 3 times (to work with ntpos table)
    new_codon_list = []  # empty list to append to
    for i in range(0, len(codon_list)): # for each of the 3 items in the codon list
        string = codon_list[i] # pull out each of the 3 repeats
        new_codon_list.append(string[:i] + string[i].upper() + string[i + 1:]) # make the ntpos of the codon pos uppercase
    return(new_codon_list) # return a new codon list to make into a df

def translate(sequence): 
    codons = []
    amino_acids = []

    if len(sequence)%3 == 0: # confirm that it can be translated
        for i in range(0, len(sequence), 3): 
            codon = sequence[i:i + 3] # pull out the codon
            codons.extend(adjust_codon(codon)) # adjust the codon to be multiples of 3 w/nts uppercase at each position
            a = aminoacid[codon]  # translate
            amino_acids.extend([a for i in range(3)])  # extend the amino acid list to have the aa listed 3 times for the df

    return(codons, amino_acids)

def adjust_AA(codon, original_aa):
    NTlist = ['A','T','G','C']
    upperIndex = getindex(codon)[0] #identify where nucleotide is
    fullcodon = list(codon) #make codon a list of three nucleotides
    currentNt = fullcodon[upperIndex]
    nonsyn_count = []
    syn_count = []

    for nt in NTlist:
        if nt != currentNt:
            adjust_codon = fullcodon
            adjust_codon[upperIndex] = nt.upper() #put minor nt in index position and make uppercase
            adjust_codon = ''.join(adjust_codon) #join into string (no longer list)
            adjustaa = aminoacid[adjust_codon.lower()]

            if adjustaa != original_aa:
                nonsyn_count.append(1/3)

            elif adjustaa == original_aa:
                syn_count.append(1/3)

    return sum(nonsyn_count), sum(syn_count)

feat = pd.read_csv(args.features,keep_default_na=False) # feature file with info on the gene must have NAME, START, END, SEGMENT
fasta = open_fasta(args.ref)

data = {
    "segment": [],
    "cds": [],
    "ntpos": [],
    "nt": [],
    "codon": [],
    "aa": []
}

for index, row in feat.iterrows():
    # pull out general information to use
    start = row['START']
    end = row['END']
    seq = fasta[row['SEGMENT']][row['START'] - 1:row['END']]    
    c, a = translate(seq.lower())  # translate

    data["segment"].extend([row['SEGMENT']] * len(seq))
    data["cds"].extend([row['NAME']] * len(seq))
    data["ntpos"].extend(range(start, end + 1))
    data["nt"].extend(seq)
    data["codon"].extend(c)
    data["aa"].extend(a)

    
df = pd.DataFrame(data)

for index, row in df.iterrows():
    n, s = adjust_AA(row['codon'], row['aa']) # goes through each position and tests 
    df['n'] = n
    df['s'] = s


print('nonsyn:')
print(df.groupby('segment')['n'].sum())

print("")
print("syn")
print(df.groupby('segment')['s'].sum())

print("")
TotalNonsyn = round(df['n'].sum())
TotalSyn = round(df['s'].sum())
Total = TotalNonsyn + TotalSyn

print("Total nonsyn: ", TotalNonsyn)
print("Total syn: ", TotalSyn)
print("Total positions (should be same as length): ", Total)
print(df.shape)

