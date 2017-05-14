#R code
#PLOT SILENT VS MISSENSE VS NONSENSE
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
require(cowplot)

flooring <- function(v){
  if (v < 0.01){return (0.01)}
  else{return (v)}
  }

SilTable <- read_tsv('result/HK68_Tlib.count') %>% 
              filter(grepl('WT-',mut)) %>%
              filter(!grepl('WT-SIL',mut)) %>%
              filter(!grepl('WT-GGTCTGAGT',mut)) %>%
              filter(InputCount>50) %>%
              .$R1Fitness %>%
              cbind(rep('SIL',length(.))) %>%
              data.frame(.) %>%
              setNames(c('Fit','Type')) %>%
              mutate(Fit=as.numeric(as.character(Fit)))
FitTable <- read_tsv('result/HK68_MutFitTable.tsv')
NonTable <- FitTable %>% 
              filter(grepl('_',Mut)) %>%
              .$Fitness %>%
              cbind(rep('NON',length(.))) %>%
              data.frame(.) %>%
              setNames(c('Fit','Type')) %>%
              mutate(Fit=as.numeric(as.character(Fit)))
MisTable <- FitTable %>%
              filter(!grepl('_',Mut)) %>%
              .$Fitness %>%
              cbind(rep('MIS',length(.))) %>%
              data.frame(.) %>%
              setNames(c('Fit','Type')) %>%
              mutate(Fit=as.numeric(as.character(Fit)))

ComTable <- rbind(MisTable,NonTable,SilTable) %>%
              mutate(Type=factor(Type,level=c('MIS','NON','SIL'))) %>%
              mutate(Fit=as.numeric(as.character(Fit))) %>%
              mutate(Fit=mapply(flooring,Fit)) %>%
              mutate(Fit=log10(Fit))

p <- ggplot(ComTable, aes(x=Type, y=Fit)) +
       geom_jitter(position=position_jitter(0.3),size=0.3,color='#56B4E9') +
       geom_violin(fill=NA,adjust=5,color='black',size=0.8) +
       geom_hline(yintercept=0,linetype=2,color='black')

ggsave('graph/HK68_FitQC_stripchart.png', p, height=3, width=2.5)
