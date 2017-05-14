#!/usr/bin/python
import os
import sys
import string
import operator
from Bio import SeqIO
from itertools import imap
from collections import Counter

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

def ID2mut(ID,WT): 
    mut   = []
    WT_resi = WT.rsplit('-')
    WT_pos  = [resi[0:-1] for resi in WT_resi]
    WT_ID   = ''.join([resi[-1] for resi in WT_resi])
    if ID == WT_ID: return 'WT'
    for pos, aa in zip(WT_pos,ID):
      resi = str(pos)+aa
      if resi in WT_resi: continue
      else: mut.append(resi)
    mut = '-'.join(mut)
    return mut

def CutoffMutFitTable(infile, outfile, WT):
  infile    = open(infile,'r')
  outfile   = open(outfile,'w')
  HK68fithash  = {}
  #outfile.write('Mut'+"\t"+'ID'+"\t"+'Fitness'+"\n")
  outfile.write('Mut'+"\t"+'Fitness'+"\n")
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    ID  = line[0]
    inputcount = int(line[2])
    Fitness  = line[7]
    if '-' in ID: continue
    if inputcount < 180: continue
    ID = 'GLS' if ID=='WT' else ID
    HK68fithash[ID]=Fitness
    mut = ID2mut(ID, WT)
    #outfile.write(mut+"\t"+ID+"\t"+Fitness+"\n") 
    outfile.write(mut+"\t"+Fitness+"\n") 
  infile.close()
  outfile.close()
  return HK68fithash

def hashin(infile,WT):
  infile = open(infile,'r')
  fithash = {}
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    fithash[mut2ID(line[0],WT)] = line[1]
  infile.close()
  return fithash

def main():
  WT = '225G-226L-228S'
  outfile = open('result/MutCompareTable.tsv','w')
  HK68fithash = CutoffMutFitTable('result/HK68_Tlib.count', 'result/HK68_MutFitTable.tsv',WT)
  WSNfithash  = hashin('result/WSN_MutFitTable.tsv','225D-226Q-228G')
  muts = set(HK68fithash.keys()).intersection(set(WSNfithash.keys()))
  print "Total number of mut in WSN: %i" % len(WSNfithash)
  print "Total number of mut in HK68: %i" % len(HK68fithash)
  print "Total number of overlapping muts: %i" % len(muts)
  outfile.write("\t".join(['mut','WSNfit','HK68fit'])+"\n")
  for mut in muts: 
    outfile.write("\t".join([mut, WSNfithash[mut], HK68fithash[mut]])+"\n")
  outfile.close()
  
if __name__ == "__main__":
  main()
