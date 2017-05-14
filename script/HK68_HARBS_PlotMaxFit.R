#R code
library(ggplot2)
library(reshape)
library(scales)
library(dplyr)
require(cowplot)

flooring <- function(v){
  if (v < 0.01){return (0.01)}
  else{return (v)}
  }

plotheatmap <- function(filename, graphname, n, h, w, cap){
  mutfit <- read.table(filename, header=1, check.names=FALSE)
  mutfit$Pos <- as.character(mutfit$Pos)
  y <- melt(mutfit)
  y$value[y$value==-1] <- NA
  y <- y[which(y$Pos>=n),]
  m <- y
  y$value[y$value>cap] <- cap
  p <- ggplot(data=y, aes(x=variable, y=Pos, fill=value)) +
              geom_tile() + 
              scale_fill_gradientn(colours=c("white", "yellow", "red"),
                         values=rescale(c(0, 0.5, 1)),
                         guide="colorbar",
                         na.value="grey") +
              theme_classic() +
              theme(panel.border = element_rect(colour = "black", fill=NA, size=4),
                    text = element_text(size=20))
  ggsave(graphname,p, height=h, width=w)
  return (m)
  }

m1 <- plotheatmap('result/HK68_MaxFitMut_1.tsv','graph/HK68_MaxFitMut_1.png',225, 2, 7, 1)
m2 <- plotheatmap('result/HK68_MaxFitMut_2.tsv','graph/HK68_MaxFitMut_2.png',225, 2, 7, 1)
m3 <- plotheatmap('result/HK68_MaxFitMut_3.tsv','graph/HK68_MaxFitMut_3.png',225, 2, 7, 1)
t1 <- m1 %>% 
        arrange(value) %>% 
        mutate(mut=mapply(function(x,y){return(paste(x,y,sep=''))},.$Pos,.$variable)) %>%
        rename(value1=value) %>%
        select(mut, value1)
t2 <- m2 %>% 
        arrange(value) %>% 
        mutate(mut=mapply(function(x,y){return(paste(x,y,sep=''))},.$Pos,.$variable)) %>%
        rename(value2=value) %>%
        select(mut, value2)
t3 <- m3 %>%
        arrange(value) %>%
        mutate(mut=mapply(function(x,y){return(paste(x,y,sep=''))},.$Pos,.$variable)) %>%
        rename(value3=value) %>%
        select(mut, value3)
t  <- t1 %>%
        inner_join(t2) %>%
        inner_join(t3) %>%
        filter(mut!='225G') %>%
        filter(mut!='226L') %>%
        filter(mut!='228S') %>%
        melt() %>%
        filter(value != 'NA') %>%
        mutate(value=mapply(flooring,.$value))
p <- ggplot(t,aes(x=variable, y=log10(value),fill=variable,color=variable)) +
       geom_jitter(width=0.2, size=0.3) +
       geom_violin(alpha=0.2)
ggsave('graph/HK68_MaxFitMut_Dist.png',p,width=3.5,height=3)
