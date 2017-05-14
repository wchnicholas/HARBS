#!/usr/bin/python
import os
import sys
import string
import operator
from Bio import SeqIO
from itertools import imap
from collections import Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def ProcessTlib(R1file):
  R2file = R1file.replace('_R1_','_R2_')
  R1RandPos = 48
  R2RandPos = 48
  Length  = 12
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  muts = []
  for R1record in R1records:
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi  = R1seq[R1RandPos:R1RandPos+Length]
    R2roi  = R2seq[R2RandPos:R2RandPos+Length]
    R1flank5 = R1seq[R1RandPos-3:R1RandPos]
    R1flank3 = R1seq[R1RandPos+Length:R1RandPos+Length+3]
    R2flank5 = R2seq[R2RandPos-3:R2RandPos]
    R2flank3 = R2seq[R2RandPos+Length:R2RandPos+Length+3]
    if R1flank5 != 'AGG' or R1flank3 != 'AGA' or R2flank5 != 'TCT' or R2flank3 != 'CCT': continue
    if 'N' in R1roi or 'N' in R2roi: continue
    if R1roi != rc(R2roi): continue
    roi  = R1roi[0:6]+R1roi[9:12]
    Tmut = translation(roi)
    if Tmut == 'GLS':
      muts.append('WT')
      muts.append('WT-'+roi)
      if roi != 'GGTCTGAGT':
        muts.append('WT-SIL')
    else: 
      muts.append(Tmut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def ProcessNA(R1file,NAoutfile,timepoint):
  R2file  = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  resi151   = []
  for R1record in R1records:
    R2record = R2records.next()
    R1_resi150 = R1record.seq[49:52]
    R1_resi151 = R1record.seq[52:55]
    R1_resi152 = R1record.seq[55:58]
    R2_resi150 = rc(R2record.seq[55:58])
    R2_resi151 = rc(R2record.seq[52:55])
    R2_resi152 = rc(R2record.seq[49:52])
    if R1_resi150==R2_resi150=='CAT' and R1_resi151==R2_resi151 and R1_resi152==R2_resi152=='AGA':
      resi151.append(translation(R1_resi151))
  resi151 = Counter(resi151)
  R1records.close()
  R2records.close()
  [NAoutfile.write(timepoint+"\t"+aa+"\t"+str(resi151[aa])+"\n") for aa in sorted(resi151.keys())]

def Output(InputDict, R1Dict, R2Dict, outfile):
  outfile = open(outfile,'w')
  muts = list(set(InputDict.keys()+R1Dict.keys()+R2Dict.keys()))
  WT_R1EnrichRatio = float(R1Dict['WT'])/float(InputDict['WT']+1)
  WT_R2EnrichRatio = float(R2Dict['WT'])/float(InputDict['WT']+1)
  outfile.write("\t".join(['mut','mutclass','InputCount','R1Count','R2Count','R1EnrichRatio','R2EnrichRatio','R1Fitness','R2Fitness'])+"\n")
  for mut in muts:
    R1EnrichRatio = float(R1Dict[mut])/float(InputDict[mut]+1)
    R2EnrichRatio = float(R2Dict[mut])/float(InputDict[mut]+1)
    R1Fitness     = float(R1EnrichRatio)/float(WT_R1EnrichRatio)
    R2Fitness     = float(R2EnrichRatio)/float(WT_R2EnrichRatio)
    mutclass    = 0 if 'WT' == mut or 'WT-' in mut else hamming('GLS',mut)
    outfile.write("\t".join(map(str,[mut,mutclass,InputDict[mut],R1Dict[mut],R2Dict[mut],R1EnrichRatio,R2EnrichRatio,R1Fitness,R2Fitness]))+"\n")
  outfile.close()

def main():
  outfile   = 'result/HK68_Tlib.count'
  InputDict = ProcessTlib('fastq/HK68-Tlib-1_S1_L001_R1_001.fastq')
  print 'Finished processing input library'
  R1Dict    = ProcessTlib('fastq/HK68-Tlib-2_S2_L001_R1_001.fastq')
  print 'Finished processing round 1 library'
  R2Dict    = ProcessTlib('fastq/HK68-Tlib-3_S3_L001_R1_001.fastq')
  print 'Finished processing round 2 library'
  Output(InputDict, R1Dict, R2Dict, outfile)

if __name__ == "__main__":
  main()
