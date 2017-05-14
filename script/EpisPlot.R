#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
require(cowplot)

freqplot <- function(StrainTable,graphname){
  p <- ggplot(StrainTable, aes(PairID,EpisCount/All)) + 
         geom_bar(stat='identity') +
         ylim(0,0.2) +
         geom_hline(yintercept=0)
  ggsave(graphname,p,height=2,width=1.5)
  }

FisherTest <- function(n1,k1,n2,k2){
  print(fisher.test(rbind(c(k1,n1-k1), c(k2,n2-k2)), alternative="two.sided"))
  }

EpisTable <- read_tsv('result/EpisCountAroundWT.tsv')
HK68 <- filter(EpisTable,Strain=='HK68')
WSN  <- filter(EpisTable,Strain=='WSN')
freqplot(HK68,'graph/HK68_EpisCount.png')
freqplot(WSN, 'graph/WSN_EpisCount.png')

HK68_225vs226 <- filter(EpisTable,Strain=='HK68',PairID=='225-226')
HK68_226vs228 <- filter(EpisTable,Strain=='HK68',PairID=='226-228')
HK68_225vs228 <- filter(EpisTable,Strain=='HK68',PairID=='225-228')
WSN_225vs226 <- filter(EpisTable,Strain=='WSN',PairID=='225-226')
WSN_226vs228 <- filter(EpisTable,Strain=='WSN',PairID=='226-228')
WSN_225vs228 <- filter(EpisTable,Strain=='WSN',PairID=='225-228')
print ('For HK68 around WT 225-226 vs 226-228')
FisherTest(HK68_225vs226$All,HK68_225vs226$EpisCount,HK68_226vs228$All,HK68_226vs228$EpisCount)
print ('For HK68 around WT 225-226 vs 225-228')
FisherTest(HK68_225vs226$All,HK68_225vs226$EpisCount,HK68_225vs228$All,HK68_225vs228$EpisCount)
print ('For HK68 around WT 225-228 vs 226-228')
FisherTest(HK68_225vs228$All,HK68_225vs228$EpisCount,HK68_226vs228$All,HK68_226vs228$EpisCount)
print ('For WSN around WT 225-226 vs 226-228')
FisherTest(WSN_225vs226$All,WSN_225vs226$EpisCount,WSN_226vs228$All,WSN_226vs228$EpisCount)
print ('For WSN around WT 225-226 vs 225-228')
FisherTest(WSN_225vs226$All,WSN_225vs226$EpisCount,WSN_225vs228$All,WSN_225vs228$EpisCount)
print ('For WSN around WT 225-228 vs 226-228')
FisherTest(WSN_225vs228$All,WSN_225vs228$EpisCount,WSN_226vs228$All,WSN_226vs228$EpisCount)
