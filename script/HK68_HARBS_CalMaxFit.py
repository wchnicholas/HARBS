#!/usr/bin/python
import os
import sys
import glob

def mutmaxfit(FitDict, viablefit, MutDict, n):
  for muts in FitDict.keys():
    if muts == 'WT': continue
    fit = FitDict[muts]
    mutcount = muts.count('-')+1
    if mutcount > n: continue
    for mut in muts.rsplit('-'):
      pos = mut[0:-1]
      aa  = mut[-1]
      if fit > MutDict[pos][aa]: MutDict[pos][aa] = fit
  return MutDict

def hashinfittable(infile):
  FitDict = {}
  infile = open(infile,'r')
  for line in infile:
    if 'Fitness' in line: continue
    line = line.rstrip().rsplit("\t")
    FitDict[line[0]] = float(line[1])
  infile.close()
  return FitDict

def genAAmat(allpos, aas):
  MutDict = {}
  for pos in allpos:
    MutDict[pos] = {}
    for aa in aas:
      MutDict[pos][aa] = -1
  return MutDict

def PrintMutDict(MutDict, aas, outfile):
  outfile = open(outfile,'w')
  allpos = sorted(MutDict.keys())
  outfile.write("\t".join(['Pos']+aas)+"\n")
  for pos in allpos:
    outfits = [pos]
    for aa in aas:
      outfits.append(MutDict[pos][aa])
    outfile.write("\t".join(map(str,outfits))+"\n")
  outfile.close()

def main():
  aas    = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_']
  WTs    = ['225G','226L','228S']
  filename  = 'result/HK68_MutFitTable.tsv'
  viablefit = 0.4
  FitDict = hashinfittable(filename)
  allpos  = sorted(list(set([mut[0:-1] for muts in FitDict.keys() for mut in muts.rsplit('-') if mut != 'WT'])))
  MutDict = genAAmat(allpos, aas)
  for WT in WTs: MutDict[WT[0:-1]][WT[-1]] = float(1)
  for n in [1,2,3]:
    outfile = 'result/HK68_MaxFitMut_'+str(n)+'.tsv'
    print "Writing %s" % outfile
    MutDict = mutmaxfit(FitDict, viablefit, MutDict, n)
    PrintMutDict(MutDict, aas, outfile)

if __name__ == "__main__":
  main()
