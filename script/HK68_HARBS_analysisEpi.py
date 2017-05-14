#!/usr/bin/python
import os
import sys
import glob

def FindEpistasis(fithash):
  muts  = fithash.keys()
  Tmuts = filter(lambda x:x.count('-')==2, muts)
  Dmuts = filter(lambda x:x.count('-')==1, muts)
  for Tmut in Tmuts: 
    mutA = Tmut.rsplit('-')[0]
    mutB = Tmut.rsplit('-')[1]
    mutC = Tmut.rsplit('-')[2]
    mutAB = mutA+'-'+mutB
    mutBC = mutB+'-'+mutC
    mutAC = mutA+'-'+mutC
    TFit  = fithash[Tmut]
    if TFit > 1:
      Interest = 0
      if mutA in muts and mutB in muts and mutAB in muts: 
        if fithash[mutA] < 0.1 and fithash[mutB] < 0.1 and fithash[mutAB] < 0.1:
          Interest = 1
      if mutB in muts and mutC in muts and mutBC in muts: 
        if fithash[mutB] < 0.1 and fithash[mutC] < 0.1 and fithash[mutBC] < 0.1:
          Interest = 1
      if mutA in muts and mutC in muts and mutAC in muts: 
        if fithash[mutA] < 0.1 and fithash[mutC] < 0.1 and fithash[mutAC] < 0.1:
          Interest = 1
      if Interest == 1:
        print Tmut, TFit
        #print "\n".join([m+"\t"+str(fithash[m]) for m in [mutA, mutB, mutC, mutAB, mutAC, mutBC] if m in fithash.keys()])
  for Dmut in Dmuts:
    mutA = Dmut.rsplit('-')[0]
    mutB = Dmut.rsplit('-')[1]
    DFit = fithash[Dmut]
    if DFit > 0.4: 
      if mutA in muts and mutB in muts:
        FitA = fithash[mutA] 
        FitB = fithash[mutB] 
        if FitA < 0.2 or FitB < 0.2: 
          print Dmut, FitA, FitB, DFit
      
def hashin(infile):
  infile = open(infile,'r')
  fithash = {}
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line = line.rstrip().rsplit("\t")
    fithash[line[0]] = float(line[1])
  infile.close()
  return fithash

def main():
  fitfile = '../../WSN/Mapping/result/MutFitTable'
  fithash = hashin(fitfile)
  FindEpistasis(fithash)
  
if __name__ == "__main__":
  main()
