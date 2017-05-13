#!/usr/bin/python
import os
import sys
import string
from Bio import SeqIO
from collections import Counter

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
  R1RandPos = 154
  R2RandPos = 114
  Length  = 18
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  muts = []
  for R1record in R1records:
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi  = R1seq[R1RandPos:R1RandPos+Length]
    R2roi  = R2seq[R2RandPos:R2RandPos+Length]
    Tmut = ''
    if 'N' in R1roi or 'N' in R2roi: continue
    if R1roi != rc(R2roi): continue
    if R1roi[0:3] == 'AGA' and R1roi[9:12] == 'CAT' and R1roi[15:18] == 'AGG': Tmut = 'WT'
    if R1roi[0:3] == 'AGG' and R1roi[9:12] == 'CAC' and R1roi[15:18] == 'AGA': Tmut = 'Mut'
    if Tmut == '': continue
    if Tmut == 'Mut':
      mut = []
      if translation(R1roi[3:6]) != 'D': mut.append('225'+translation(R1roi[3:6])) 
      if translation(R1roi[6:9]) != 'Q': mut.append('226'+translation(R1roi[6:9])) 
      if translation(R1roi[12:15]) != 'G': mut.append('228'+translation(R1roi[12:15])) 
      if len(mut) == 0:
        mut = 'SIL'
      else:
        mut = '-'.join(mut)
    if Tmut == 'WT':
      mut = 'WT'
    muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def ProcessTlib2(R1file):
  R2file = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  seglen    = 200
  muts = []
  for R1record in R1records:
    R2record  = R2records.next()
    R1roi  = str(R1record.seq)[0:200]
    R2roi  = str(R2record.seq)[0:200]
    Tmut = ''
    codonBC = R1roi[49:52]+R1roi[55:58]+R1roi[154:157]+R1roi[163:166]
    if 'N' in R1roi or 'N' in R2roi: continue
    if R1roi != rc(R2roi): continue
    if codonBC == 'GACCAGAGGCAC': Tmut = 'Mut'
    if codonBC == 'GATCAAAGACAT': Tmut = 'WT'
    if Tmut == '': continue
    if Tmut == 'Mut':
      mut = []
      if translation(R1roi[52:55]) != 'E': mut.append('190'+translation(R1roi[52:55])) 
      if translation(R1roi[157:160]) != 'D': mut.append('225'+translation(R1roi[157:160])) 
      if translation(R1roi[160:163]) != 'Q': mut.append('226'+translation(R1roi[160:163])) 
      if len(mut) == 0:
        mut = 'SIL'
      else:
        mut = '-'.join(mut)
    if Tmut == 'WT':
      mut = 'WT'
    muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def ProcessDlib(R1file):
  R2file = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file, "fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  muts = []
  for R1record in R1records:
    R2record = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1resi134 = R1seq[33:42]
    R1resi136 = R1seq[39:48]
    R1resi153 = R1seq[90:99];   R2resi153 = rc(R2seq[253:262])
    R1resi155 = R1seq[96:105];  R2resi155 = rc(R2seq[247:256])
    R1resi183 = R1seq[180:189]; R2resi183 = rc(R2seq[163:172])
    R1resi190 = R1seq[201:210]; R2resi190 = rc(R2seq[142:151])
    R1resi194 = R1seq[213:222]; R2resi194 = rc(R2seq[130:139])
    R1resi195 = R1seq[216:225]; R2resi195 = rc(R2seq[127:136])
    R2resi225 = rc(R2seq[37:46])
    R2resi226 = rc(R2seq[34:43])
    R2resi228 = rc(R2seq[28:37])
    mutresi = []
    if R1resi134[0:3]=='AAT' and R1resi134[6:9]=='GTT' and 'N' not in R1resi134:
      mutresi.append('134'+translation(R1resi134[3:6]))
      if R1resi136[6:9]=='GTC' and 'N' not in R1resi136: mutresi.append('136'+translation(R1resi136[3:6])) 
    if R1resi136[0:3]=='GTC' and R1resi136[6:9]=='GTC' and 'N' not in R1resi136:
      mutresi.append('136'+translation(R1resi136[3:6]))
    if R1resi153[0:3]=='CTT' and R1resi153[6:9]=='CTC' and R1resi153==R2resi153 and 'N' not in R1resi153:
      mutresi.append('153'+translation(R1resi153[3:6]))
      if R1resi155[6:9]=='AAA' and 'N' not in R1resi155: mutresi.append('155'+translation(R1resi155[3:6]))
    if R1resi155[0:3]=='CTA' and R1resi155[6:9]=='AAA' and R1resi155==R2resi155 and 'N' not in R1resi155:
      mutresi.append('155'+translation(R1resi155[3:6]))
    if R1resi183[0:3]=='GTA' and R1resi183[6:9]=='CAT' and R1resi183==R2resi183 and 'N' not in R1resi183:
      mutresi.append('183'+translation(R1resi183[3:6]))
    if R1resi190[0:3]=='GAC' and R1resi190[6:9]=='CAG' and R1resi190==R2resi190 and 'N' not in R1resi190:
      mutresi.append('190'+translation(R1resi190[3:6]))
    if R1resi194[0:3]=='AGC' and R1resi194[6:9]=='TAC' and R1resi194==R2resi194 and 'N' not in R1resi194:
      mutresi.append('194'+translation(R1resi194[3:6]))
    if R1resi195[0:3]=='CTA' and R1resi195[6:9]=='AGC' and R1resi195==R2resi195 and 'N' not in R1resi195:
      mutresi.append('195'+translation(R1resi195[3:6]))
    if R2resi225[0:3]=='AGG' and R2resi225[6:9]=='CAG' and 'N' not in R2resi225:
      mutresi.append('225'+translation(R2resi225[3:6]))
    if R2resi226[0:3]=='GAC' and R2resi226[6:9]=='CAC' and 'N' not in R2resi226:
      mutresi.append('226'+translation(R2resi226[3:6]))
    if R2resi228[0:3]=='CAC' and R2resi228[6:9]=='AGA' and 'N' not in R2resi228:
      mutresi.append('228'+translation(R2resi228[3:6]))
      if R2resi225[0:3]=='AGG' and R2resi225[6:9]=='CAA' and 'N' not in R2resi225:
        mutresi.append('225'+translation(R2resi225[3:6]))
    if R1resi194[0:3]=='AGC' and R1resi195[6:9]=='AGC' and 'N' not in R1resi194 and 'N' not in R1resi195:
      mutresi.append('194'+translation(R1resi194[3:6]))
      mutresi.append('195'+translation(R1resi195[3:6]))
    if R2resi225[0:3]=='AGG' and R2resi226[6:9]=='CAC' and 'N' not in R2resi225 and 'N' not in R2resi226:
      mutresi.append('225'+translation(R2resi225[3:6]))
      mutresi.append('226'+translation(R2resi226[3:6]))
    if len(mutresi) > 2: continue
    if len(mutresi) == 0: 
      if ( R1resi134[0:3]=='AAC' and R1resi134[6:9]=='GTA' and R1resi136[6:9]=='GTA' and R1resi136[6:9]=='GTA' and 
           R1resi153[0:3]=='CTA' and R1resi153[6:9]=='CTG' and R1resi155[0:3]=='CTG' and R1resi155[6:9]=='AAG' and 
           R1resi183[0:3]=='GTT' and R1resi183[6:9]=='CAC' and R1resi190[0:3]=='GAT' and R1resi190[6:9]=='CAA' and 
           R1resi194[0:3]=='AGT' and R1resi194[6:9]=='TAT' and R1resi195[0:3]=='CTC' and R1resi195[6:9]=='AGT' and 
           R2resi225[0:3]=='AGA' and R2resi225[6:9]=='CAA' and R2resi226[0:3]=='GAT' and R2resi226[6:9]=='CAT' and 
           R2resi228[0:3]=='CAT' and R2resi228[6:9]=='AGG' ):
        mut = 'WT'
        muts.append(mut)
    else: 
      mutresiadj=list(set(mutresi).difference(set(['134G','136T','153W','155T','183H',
                                                   '190E','194L','195Y','225D','226Q','228G'])))
      if len(mutresiadj) > 2: continue
      elif len(mutresiadj) == 0: mut = 'SIL'+':'+"-".join(mutresi)
      elif len(mutresiadj) == 1: mut = mutresiadj[0]
      elif len(mutresiadj) == 2: mut = "-".join(sorted(mutresiadj,key=lambda x:int(x[0:3])))
      else: 'Something is wrong!'; sys.exit()
      muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def Output(InputDict, R2Dict, outfile):
  outfile = open(outfile,'w')
  muts = list(set(InputDict.keys()+R2Dict.keys()))
  WTEnrichRatio = float(R2Dict['WT'])/float(InputDict['WT']+1)
  outfile.write("\t".join(['mut','mutclass','InputCount','R2Count','EnrichRatio','Fitness'])+"\n")
  for mut in muts:
    EnrichRatio = float(R2Dict[mut])/float(InputDict[mut]+1)
    Fitness     = float(EnrichRatio)/float(WTEnrichRatio)
    mutclass    = 0 if 'WT' in mut or 'SIL' in mut else mut.count('-')+1
    outfile.write("\t".join(map(str,[mut,mutclass,InputDict[mut],R2Dict[mut],EnrichRatio,Fitness]))+"\n")
  outfile.close()

def wrapper(libinputfile, libR2file, outfile, ProcessFunction):
  InputDict = ProcessFunction(libinputfile)
  print 'Finish processing %s' % libinputfile
  R2Dict    = ProcessFunction(libR2file)
  print 'Finish processing %s' % libR2file
  Output(InputDict, R2Dict, outfile)
  print 'Finish writing %s' % outfile

def main():
  wrapper('fastq/WSN_HARBS-1_S12_L001_R1_001.fastq','fastq/WSN_HARBS-4_S15_L001_R1_001.fastq','result/WSN_SingleMutLib.count',ProcessDlib)
  wrapper('fastq/HARBS-5_S16_L001_R1_001.fastq','fastq/WSN_HARBS-2_S13_L001_R1_001.fastq','result/WSN_DoubleMutLib.count',ProcessDlib)
  wrapper('fastq/WSN_HARBS-3_S14_L001_R1_001.fastq','fastq/WSN_HARBS-6_S17_L001_R1_001.fastq','result/WSN_TripleMutLib.count',ProcessTlib)

if __name__ == "__main__":
  main()
