#!/usr/bin/env python3
import os
import math
import argparse
import pandas as pd
import sys

#    merge_freqandbim.py  --freq  ${headout}_tmp.frq --bim $bim --out ${headout}.frq
def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--out', type=str,help="out of tex file",default="test.tex")
    parser.add_argument('--freq',type=str,required=False,help="beta header in inp files")
    parser.add_argument('--bim',type=str,required=False,help="beta header in inp files")
    args = parser.parse_args()
    return args

args = parseArguments()

# CHR               SNP   A1   A2          MAF  NCHROBS
freq= pd.read_csv(args.freq,delim_whitespace=True)
#1	1:723918:G:A	0.417725	723918	A	G

bim=pd.read_csv(args.bim,delim_whitespace=True,header=None)
bim.columns= ["CHR", "SNP", "CM", "BP", "A1","A2"]
bim=bim[['CHR', 'SNP', 'BP']]

allfreq=freq.merge(bim, left_on=['CHR','SNP'], right_on=['CHR','SNP'])
allfreq.to_csv(args.out,sep="\t",header=True,index=False,na_rep="NA")




