#!/usr/bin/python
import os
import sys
import glob
import operator
from itertools import imap
from collections import defaultdict

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def mut2ID(mut, WT):
  WT_resi = WT.rsplit('-')
  if mut == 'WT': 
    return ''.join([resi[-1] for resi in WT_resi])
  for resi in WT_resi:
    pos = resi[0:-1]
    aa  = resi[-1]
    if pos not in mut:
      mut += '-'+resi
  mut = sorted(mut.rsplit('-'),key=lambda x:x[0:-1])
  ID  = ''.join([resi[-1] for resi in mut])
  return ID

def fithashing(fitfile, WT):
  infile  = open(fitfile,'r')
  fitdict = {}
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = float(line[1])
    fitdict[mut] = fit
  infile.close()
  return fitdict

def fitclean(fitdict, WT):
  WTpos = set([resi[0:-1] for resi in WT.rsplit('-')])
  for mut in fitdict.keys():
    fit = fitdict[mut]
    del fitdict[mut]
    if '_' in mut: continue
    pos = set([resi[0:-1] for resi in mut.rsplit('-')])
    if pos.issubset(WTpos) or mut == 'WT': 
      ID  = mut2ID(mut,WT)
      fitdict[ID] = fit
  return fitdict

def ProcessPair(mut1, mut2):
  mut3 = ''
  mut4 = ''
  match  = ''
  switch = 0
  for n in range(len(mut1)):
    if mut1[n] == mut2[n]:
      match = n
      mut3 += mut1[n]
      mut4 += mut1[n]
    elif switch == 0: 
      mut3 += mut1[n]
      mut4 += mut2[n]
      switch += 1
    elif switch == 1:
      mut3 += mut2[n]
      mut4 += mut1[n]
  return mut3, mut4, match

def writeEpiMap(Epi_AllPairDict, Epi_AAPairDict, outfile, WT):
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  resis  = [pos[-1]+pos[0:-1]+aa for pos in WT.rsplit('-') for aa in aas if pos[-1] != aa]
  outfile = open(outfile,'w')
  header = "\t".join(['mut1','mut2','value'])
  outfile.write(header+"\n")
  for resi1 in resis:
    for resi2 in resis:
      ID = '-'.join(sorted([resi1,resi2],key=lambda x:int(x[1:-1])))
      if Epi_AllPairDict[ID]==0: value = -1
      elif Epi_AllPairDict[ID]>0 and Epi_AAPairDict[ID]==0: value = 0
      elif Epi_AllPairDict[ID]>0 and Epi_AAPairDict[ID]>0:  value = float(Epi_AAPairDict[ID]/Epi_AllPairDict[ID])
      else: print "Something is wrong with Epi_AAPair classification"; sys.exit()
      outfile.write("\t".join([resi1, resi2,str(value)])+"\n")
  outfile.close()

def CountEpis(fitdict, WT, outfile, epifile, Strain):
  muts  = [mut for mut in fitdict.keys() if '_' not in mut]
  WT_ID = mut2ID(WT, WT)
  TotalCount = 0
  AllPairDict   = defaultdict(int)
  EpiPairDict   = defaultdict(int)
  Epi_AAPairDict    = defaultdict(int)
  All_AAPairDict    = defaultdict(int)
  print "Processing a total of %i mutants" % len(muts)
  for i in range(len(muts)): 
    for j in range(len(muts)):
      mut1 = muts[i]
      mut2 = muts[j]
      #if i > j: continue
      if mut1 != WT_ID: continue
      if hamming(mut1,mut2) == 2:
        mut3, mut4, match = ProcessPair(mut1,mut2)
        if match == 0:   pair_ID = '226-228'; pair_aa = mut1[1]+'226'+mut2[1]+'-'+mut1[2]+'228'+mut2[2]
        elif match == 1: pair_ID = '225-228'; pair_aa = mut1[0]+'225'+mut2[0]+'-'+mut1[2]+'228'+mut2[2]
        elif match == 2: pair_ID = '225-226'; pair_aa = mut1[0]+'225'+mut2[0]+'-'+mut1[1]+'226'+mut2[1]
        else: print 'Something is wrong with pair_ID'; sys.exit()
        if mut1 in muts and mut2 in muts and mut3 in muts and mut4 in muts:
          TotalCount += 1
          AllPairDict[pair_ID] += 1
          All_AAPairDict[pair_aa] += 1
          mut1fit = fitdict[mut1]
          mut2fit = fitdict[mut2]
          mut3fit = fitdict[mut3]
          mut4fit = fitdict[mut4]
          if mut3fit<mut1fit and mut4fit<mut1fit and mut3fit<mut2fit and mut4fit<mut2fit:
            EpiPairDict[pair_ID] += 1
            Epi_AAPairDict[pair_aa] += 1
  
  writeEpiMap(All_AAPairDict, Epi_AAPairDict, outfile, WT)
  
  for pair_ID in AllPairDict.keys():
     EpiPairCount = EpiPairDict[pair_ID]
     AllPairCount = AllPairDict[pair_ID]
     frac = float(EpiPairCount)/float(AllPairCount)
     print "Epistasis in position pair %s: %i out of %i (%f)" % (pair_ID, EpiPairCount, AllPairCount, frac)
     epifile.write("\t".join(map(str,[pair_ID, AllPairCount, EpiPairCount, Strain]))+"\n")
     
def wrapper(fitfile,WT,outfile,epifile,Strain):
  fitdict = fithashing(fitfile, WT)
  print "Read a total of %i mutants" % len(fitdict.keys())
  fitdict = fitclean(fitdict, WT)
  print "A total of %i clean mutants" % len(fitdict.keys())
  CountEpis(fitdict, WT, outfile, epifile, Strain)

def main():
  epifile = open('result/EpisCountAroundWT.tsv','w')
  epifile.write("\t".join(['PairID','All','EpisCount','Strain'])+"\n")
  wrapper('result/HK68_MutFitTable.tsv', '225G-226L-228S', 'result/HK68_EpisMap.tsv', epifile, 'HK68')
  wrapper('result/WSN_MutFitTable.tsv', '225D-226Q-228G', 'result/WSN_EpisMap.tsv', epifile, 'WSN')
  epifile.close()
  
  
if __name__ == "__main__":
  main()
