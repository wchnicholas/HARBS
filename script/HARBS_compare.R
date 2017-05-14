#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
require(cowplot)

floor <- function(v,fvalue){
  v <- as.numeric(v)
  if (v < fvalue){return (fvalue)}
  else{return (v)}
  }

coloring <- function(mut){
  WT <- c('DQG','GLS')
  WSNValid  <- c('LQG','EQG','QQG','GQG','DSG','DNG','DMG','DAG','DHG',
                 'DQA','DQS','HQG','EMG','EAG','ENG','QMG','QAG','QHG',
                 'LSG','EQA','LQS','DNA','DAA','DSS','LSS','FSS','ENA',
                 'EAA','MTA')
  HK68Valid <- c('GQG','MTA','LSS','SWA','WGS','RAA','QAS','ENA','GNS',
                 'GLA','GNA','ELS','GRS','GLE','ERS','ELE','GRE','ERE')
  BothValid <- intersect(WSNValid,HK68Valid)
  if (grepl('_',mut)){return ('non')}
  else if (mut %in% WT){return ('WT')}
  else if (mut %in% BothValid){return ('Both')}
  else if (mut %in% WSNValid){return ('WSN')}
  else if (mut %in% HK68Valid){return ('HK68')}
  else{return ('mis')}
  }

sizing <- function(color){
  if (color=='Both' | color=='WSN' | color=='HK68' | color=='WT'){return (1.2)}
  else{return (0.8)}
  }

labeling <- function(color,mut){
  if (color=='Both' | color=='WSN' | color=='HK68' | color=='WT'){return (mut)}
  else{return ('')}
  }

FitTable <- read_tsv('result/MutCompareTable.tsv') %>%
              mutate(HK68fit=mapply(floor,HK68fit,min(HK68fit[which(HK68fit!=0)]))) %>%
              mutate(WSNfit=mapply(floor,WSNfit,min(WSNfit[which(WSNfit!=0)]))) %>%
              mutate(color=factor(mapply(coloring,mut),levels=c('mis','non','Both','WSN','HK68','WT'))) %>%
              mutate(size=mapply(sizing, color)) %>%
              mutate(label=mapply(labeling,color,mut))
FitTable_Valid <- FitTable %>%
                    filter(size==1.2)

p <- ggplot(FitTable, aes(log10(WSNfit),log10(HK68fit),color=color,size=size,label=label)) +
       geom_point() +
       scale_color_manual(values=c('#BBBBBB7F','#FF0000AF','#06cc00AF','#0081CCAF','#A0A000AF','#7001BAAF')) + 
       scale_size(range = c(0.05, 0.8)) + 
       xlim(-5.3,0.5) +
       ylim(-3.3,1.3) +
       geom_text_repel(data=FitTable_Valid,aes(log10(WSNfit),log10(HK68fit),color=color,label=label),
         size = 1.8,
         fontface = 'bold',
         segment.size = 0.2,
         point.padding = unit(0.2, "lines"),
         force = 0.5,
         nudge_x=ifelse(log10(FitTable_Valid$WSNfit)< -5 | FitTable_Valid$mut=='EQG' ,0.5,0),
         nudge_y=ifelse(FitTable_Valid$mut=='ERE',-0.05,0),
         show.legend = FALSE)
ggsave('graph/MutFitCompare.png',p,height=2.5,width=3.5)
