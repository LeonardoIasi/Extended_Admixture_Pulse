# Icelandig Fragments #
# translate bp length into genetic length using the Icelandic map #
# Fit the model #
suppressPackageStartupMessages({
  library(VGAM)
  library(tidyverse)
  library("DEoptim")
  library("MASS")
  library(bbmle)
  library(rethinking)
})

setwd("Desktop/ATE_Paper/Admixture_Time_Inference_Paper/Real_Data_Analysis/")
Fragments <- read.table("Archaicfragments_Iceland.txt",header = T)

Recomb_Map <- read.table("RecombmapIceland.txt", header = T)

Fragments_cM <- c()
for(chr in 1:23){
  if(chr==23){
    chr="X"
  }
  else{
    chr=chr
  }
  Recomb_Map_chr <- Recomb_Map[Recomb_Map$Chr==paste("chr",chr,sep=""),]
  Fragments_chr <- Fragments[Fragments$chrom==chr,]
  Start_Arch_Freg <- approx(x=Recomb_Map_chr$Begin, y=Recomb_Map_chr$cM, xout=Fragments_chr$start,method = "constant")
  Stop_Arch_Freg <- approx(x=Recomb_Map_chr$Begin, y=Recomb_Map_chr$cM, xout=Fragments_chr$end,method = "constant")
  Chr_x <- data.frame(chr=chr,Start_cM=Start_Arch_Freg$y,End_cM=Stop_Arch_Freg$y,length_bp=Fragments_chr$length,called=Fragments_chr$called)
  Fragments_cM <- rbind(Fragments_cM,Chr_x)
}

Fragments_cM <- c()
for(chr in 1:23){
  if(chr==23){
    chr="X"
  }
  else{
    chr=chr
  }
  Recomb_Map_chr <- Recomb_Map[Recomb_Map$Chr==paste("chr",chr,sep=""),]
  Fragments_chr <- Fragments[Fragments$chrom==chr,]
  Start_Arch_Freg <- approx(x=Recomb_Map_chr$Begin, y=Recomb_Map_chr$cMperMb, xout=Fragments_chr$start,method = "constant")
  Fragment_in_cM <- 
  Chr_x <- data.frame(chr=chr,Start_cM=Start_Arch_Freg$y,End_cM=Stop_Arch_Freg$y,length_bp=Fragments_chr$length,called=Fragments_chr$called)
  Fragments_cM <- rbind(Fragments_cM,Chr_x)
}

Fragments_cM_clean <- Fragments_cM[complete.cases(Fragments_cM$Start_cM),]
Fragments_cM_clean <- Fragments_cM_clean[complete.cases(Fragments_cM_clean$End_cM),]
Fragments_cM_clean$length_cM <- Fragments_cM_clean$End_cM - Fragments_cM_clean$Start_cM

cutoff=min(Fragments_cM_clean$length_cM)


plot(Fragments_cM_clean$length_cM,Fragments_cM_clean$length_bp,pch=19,col=col.alpha("black",0.2),xlab="Fragments in cM",ylab="Fragments in bp")
barplot(Fragments_cM_clean$length_bp,Fragments_cM_clean$called)
barplot(Fragments_cM_clean$length_cM,names.arg = Fragments_cM_clean$called,las=2)

Fragment_model_1 <- ulam(
  alist(
    length_cM ~ dexp( lambda) ,
    lambda <- a ,
    a ~ dexp(1)
  ) , data=Fragments_cM_clean, chains = 4,cores = 4)

precis(Fragment_model_1,prob = 0.95)

Fragment_model_2 <- ulam(
  alist(
    length_cM ~ dexp( lambda) ,
    lambda <- a ,
    a ~ dnorm(1000,500)
  ) , data=Fragments_cM_clean, chains = 4,cores = 4)

precis(Fragment_model_2,prob = 0.95)




