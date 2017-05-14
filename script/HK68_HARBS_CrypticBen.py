#!/usr/bin/python
import os
import sys
import glob
import colorsys
import networkx as nx
import operator
from itertools import imap

def hashinfittable(infile):
  FitDict = {}
  infile = open(infile,'r')
  for line in infile:
    if 'Fitness' in line: continue
    line = line.rstrip().rsplit("\t")
    FitDict[line[0]] = float(line[1])
  infile.close()
  return FitDict
  
def CrypticBenDicting(resiofinterest, aas):
  CrypticBenDict = {}
  for resi in resiofinterest:
    for aa in aas:
      CrypticBenDict[resi+aa] = 0
  return (CrypticBenDict)

def CrypticBenAnalysis(FitDict,CrypticBenDict):
  for SingleMut in CrypticBenDict.keys():
    for mut2 in FitDict.keys():
      if SingleMut in mut2:
        mut1 = mut2.rsplit('-')
        mut1.remove(SingleMut)
        if len(mut1)==0: mut1 = 'WT'
        else: mut1 = '-'.join(mut1)
        if FitDict.has_key(mut1):
          SingleMutFit = (FitDict[mut2]-FitDict[mut1])
          if CrypticBenDict[SingleMut] < SingleMutFit:
            CrypticBenDict[SingleMut] = SingleMutFit
  return (CrypticBenDict)

def CrypticBenOut(CrypticBenDict, CrypticBen_outfile, wtaa):
  outfile = open(CrypticBen_outfile,'w')
  outfile.write("\t".join(['Position','Amino Acid','CrypFitDiff'])+"\n")
  for SingleMut in CrypticBenDict.keys():
    pos = SingleMut[0:-1]  
    aa  = SingleMut[-1]
    fit = -1 if SingleMut in wtaa else CrypticBenDict[SingleMut]
    outfile.write("\t".join(map(str,[pos, aa, fit]))+"\n")
  outfile.close()

def main():
  FitDict        = hashinfittable('result/HK68_MutFitTable.tsv')
  CrypticBen_outfile = 'result/HK68_crypticben.tsv'
  viablecutoff   = 0.4
  resiofinterest = ['225','226','228']
  wtaa = ['225G','226L','228S']
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_']
  CrypticBenDict = CrypticBenDicting(resiofinterest, aas)
  CrypticBenDict = CrypticBenAnalysis(FitDict,CrypticBenDict)
  CrypticBenOut(CrypticBenDict, CrypticBen_outfile, wtaa)
  
 
if __name__ == "__main__":
  main()
