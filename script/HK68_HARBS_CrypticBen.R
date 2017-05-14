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
library(gridExtra)
library(grid)

aas <- c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_')
CrypBen_table <- read_tsv('result/HK68_crypticben.tsv') %>%
                   mutate(Position=as.character(Position)) %>%
                   mutate(`Amino Acid`=factor(`Amino Acid`,aas))
CrypBen_table$CrypFitDiff[CrypBen_table$CrypFitDiff==-1] <- NA
CrypBen_table$CrypFitDiff[CrypBen_table$CrypFitDiff>1] <- 1
p <- ggplot(data=CrypBen_table, aes(x=`Amino Acid`, y=Position, fill=CrypFitDiff)) +
            geom_tile() + 
	    scale_fill_gradientn(colours=c("white", "white", "purple"),
            values=rescale(c(0, 0.5, 1)),
		       guide="colorbar",
		       na.value="grey") +
	    theme_classic() +
	    theme(panel.border = element_rect(colour = "black", fill=NA, size=4),
		  text = element_text(size=20))
ggsave('graph/HK68_CrypticBen.png',p, height=2, width=7)
