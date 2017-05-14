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

PlotEpiMap <- function(filename,graphname,resi){
  levels <- c()
  for (pos in resi){
    for (aa in c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W')){
      levels <- c(levels,(paste(pos,aa,sep='')))
      }
    }

  t <- read_tsv(filename)
  t$value[t$value==-1] <- NA
  t <- t %>% 
         mutate(mut1=factor(mut1,levels=levels)) %>%
         mutate(mut2=factor(mut2,levels=levels))
  p <- ggplot(data=t, aes(x=mut1, y=mut2, fill=value)) +
            geom_tile(colour = "black") +
            scale_fill_gradientn(colours=c("white", "yellow", "red"),
                       values=rescale(c(0, 0.5, 1)),
                       guide='none',
                       na.value="grey") +
            theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                  text = element_text(size=6), 
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.5),
                  axis.text.y = element_text(size = 3.5)) +
            xlab('Initial Amino Acid') +
            ylab('Targeting Amino Acid')
  ggsave(graphname,p,height=3,width=3)
  }


PlotEpiMap('result/HK68_EpisMap.tsv','graph/HK68_EpiMap.png',c('G225','L226','S228'))
PlotEpiMap('result/WSN_EpisMap.tsv','graph/WSN_EpiMap.png',c('D225','Q226','G228'))
