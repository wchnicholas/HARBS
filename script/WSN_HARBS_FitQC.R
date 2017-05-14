#R code
library(ggplot2)
require(cowplot)

flooring <- function(v){
  if (v < 0.001){return (0.001)}
  else{return (v)}
  }

coloring <- function(v){
  if (v >= 0){return ('green')}
  else{return ('black')}
  }

FitnessQCplot <- function(Single, Double, Triple, graphname, h, w){
  SingleMis <- Single[which(Single$InputCount >= 5 & Single$mutclass<=1 & ! grepl('SIL',Single$mut) & ! grepl('WT',Single$mut)),]
  DoubleMis <- Double[which(Double$InputCount >= 5 & Double$mutclass<=2 & ! grepl('SIL',Double$mut) & ! grepl('WT',Double$mut)),]
  TripleMis <- Triple[which(Triple$InputCount >= 5 & Triple$mutclass<=3 & ! grepl('SIL',Triple$mut) & ! grepl('WT',Triple$mut)),]
  SingleNon <- Single[which(Single$InputCount >= 5 & Single$mutclass<=1 & grepl('_',Single$mut)),]
  DoubleNon <- Double[which(Double$InputCount >= 5 & Double$mutclass<=2 & grepl('_',Double$mut)),]
  TripleNon <- Triple[which(Triple$InputCount >= 5 & Triple$mutclass<=3 & grepl('_',Triple$mut)),]
  SingleSil <- Single[which(Single$InputCount >= 5 & Single$mutclass<=1 & grepl('SIL',Single$mut) & ! grepl('-',Single$mut)),]
  DoubleSil <- Double[which(Double$InputCount >= 5 & Double$mutclass<=2 & grepl('SIL',Double$mut)),]
  TripleSil <- Triple[which(Triple$InputCount >= 5 & Triple$mutclass<=3 & grepl('SIL',Triple$mut)),]
  SingleMis <- SingleMis$Fitness
  DoubleMis <- DoubleMis$Fitness
  TripleMis <- TripleMis$Fitness
  SingleNon <- SingleNon$Fitness
  DoubleNon <- DoubleNon$Fitness
  TripleNon <- TripleNon$Fitness
  SingleSil <- SingleSil$Fitness
  DoubleSil <- DoubleSil$Fitness
  TripleSil <- TripleSil$Fitness
  MisFit    <- log10(mapply(flooring,c(SingleMis,DoubleMis,TripleMis)))
  NonFit    <- log10(mapply(flooring,c(SingleNon,DoubleNon,TripleNon)))
  SilFit    <- log10(mapply(flooring,c(SingleSil,DoubleSil,TripleSil)))
  type      <- c(rep('Silent',length(SilFit)),rep('Nonsense',length(NonFit)),rep('Missense',length(MisFit)))
  Fit       <- c(SilFit, NonFit, MisFit)
  cat       <- c(rep('Single',length(SingleSil)),rep('Double',length(DoubleSil)),rep('Triple',length(TripleSil)),
                 rep('Single',length(SingleNon)),rep('Double',length(DoubleNon)),rep('Triple',length(TripleNon)),
                 rep('Single',length(SingleMis)),rep('Double',length(DoubleMis)),rep('Triple',length(TripleMis)))
  t <- data.frame(cat, type, Fit)
  t$cat <- factor(t$cat,levels=c('Single','Double','Triple'))
  t <- t[sample(nrow(t)),]
  p <- ggplot(t, aes(x=type, y=Fit)) + 
         geom_jitter(position=position_jitter(0.3),size=0.3, aes(color=cat)) +
         geom_violin(fill=NA,adjust=5) +
         geom_hline(yintercept=0,linetype=2,color='black')
  ggsave(graphname,p, height=h, width=w)
  }

Single <- read.table('result/WSN_SingleMutLib.count', header=1)
Double <- read.table('result/WSN_DoubleMutLib.count', header=1)
Triple <- read.table('result/WSN_TripleMutLib.count', header=1)
FitnessQCplot(Single, Double, Triple, 'graph/WSN_FitQC_stripchart.png', 3, 2)
