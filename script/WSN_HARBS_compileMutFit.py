#!/usr/bin/python
import os
import sys
import numpy as np

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

def MutHashing(infile,inputcountcutoff,outputcountcutoff,mutclasscutoff):
  mutdict = {}
  infile = open(infile,'r')
  countline = 0
  for line in infile.xreadlines():
    countline += 1
    if countline == 1: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    mutclass = int(line[1])
    fit      = float(line[5])
    inputcount  = int(line[2])
    outputcount = int(line[3])
    if mutclass > mutclasscutoff: continue
    if inputcount < inputcountcutoff: continue
    if outputcount < outputcountcutoff: continue
    mutdict[mut] = fit
  infile.close()
  return mutdict

def Compiling(SingleDict, DoubleDict, TripleDict, outfile):
  muts = list(set(SingleDict.keys()+DoubleDict.keys()+TripleDict.keys()))
  outfile = open(outfile,'w')
  #outfile.write('Mut'+"\t"+'ID'+"\t"+'Fitness'+"\n")
  outfile.write('Mut'+"\t"+'Fitness'+"\n")
  for mut in muts: 
    if 'SIL' not in mut: 
      fit = []
      if mut in SingleDict.keys(): fit.append(SingleDict[mut])
      if mut in DoubleDict.keys(): fit.append(DoubleDict[mut])
      if mut in TripleDict.keys(): fit.append(TripleDict[mut])
      ID = mut2ID(mut,'225D-226Q-228G') if mut.count('225')+mut.count('226')+mut.count('228') == mut.count('-')+1 else '-'
      #outfile.write(mut+"\t"+ID+"\t"+str(np.mean(fit))+"\n")
      outfile.write(mut+"\t"+str(np.mean(fit))+"\n")
  outfile.close()

def main():
  inputcountcutoff  = 5
  outputcountcutoff = 0
  outfile           = 'result/WSN_MutFitTable.tsv'
  SingleDict = MutHashing('result/WSN_SingleMutLib.count', inputcountcutoff, outputcountcutoff, 1)
  DoubleDict = MutHashing('result/WSN_DoubleMutLib.count', inputcountcutoff, outputcountcutoff, 2)
  TripleDict = MutHashing('result/WSN_TripleMutLib.count', inputcountcutoff, outputcountcutoff, 3)
  print 'Finish reading all count tables'
  Compiling(SingleDict, DoubleDict, TripleDict, outfile)
  

if __name__ == "__main__":
  main()
