#!/usr/bin/python
import os
import sys
import glob

def hashinfittable(infile):
  FitDict = {}
  infile = open(infile,'r')
  for line in infile:
    line = line.rstrip().rsplit("\t")
    FitDict[line[0]] = float(line[1])
  infile.close()
  return FitDict

def DoubleMutMelt(FitDict, mutIDs, outfile, WTs):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut1','mut2','fit'])+"\n")
  for mutID1 in mutIDs:
    for mutID2 in mutIDs:
      pos1 = mutID1[1::]
      pos2 = mutID2[1::]
      aa1  = mutID1[0]
      aa2  = mutID2[0]
      if pos1+aa1 in WTs: continue
      if pos2+aa2 in WTs: continue
      if int(pos1) < int(pos2): #ID IN result/WSNWT_MutFitTable ALWAYS ARRANGED BY ORDERING THE POSITIONS OF MUTATED RESIDUES
        DmutID = pos1+aa1+'-'+pos2+aa2
        fit    = FitDict[DmutID] if FitDict.has_key(DmutID) else -1
        outfile.write(pos1+aa1+"\t"+pos2+aa2+"\t"+str(fit)+"\n")
      if int(pos2) > int(pos1):
        #assert(mutID2+'-'+mutID1 not in FitDict.keys())
        DmutID = pos2+aa2+'-'+pos1+aa1
        fit    = FitDict[DmutID] if FitDict.has_key(DmutID) else -1
  outfile.close()

def genmutID(WTs):
  residues = [WT[0:-1] for WT in WTs]
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  return [aa+residue for residue in residues for aa in aas]

def main():
  filename = 'result/WSNWT_MutFitTable'
  outfile  = 'result/DmutMelt'
  WTs = '134G','136T','153W','155T','183H','190E','194L','195Y','225D','226Q','228G'
  mutIDs = genmutID(WTs)
  FitDict = hashinfittable('result/MutFitTable')
  DoubleMutMelt(FitDict, mutIDs, outfile, WTs)

if __name__ == "__main__":
  main()
