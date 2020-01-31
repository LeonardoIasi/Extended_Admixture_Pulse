### Read in MS Command and translat into MSPrime ###
args <- commandArgs(TRUE)
MS_input <- as.character(args[1])
Mut_rate <- as.numeric(args[2])
MS_Code <- paste(readLines(con = MS_input), collapse = " ")

Mut_rate_Vindjia <- 0.5e-8
library(stringr)

MS_Code <- paste(readLines(con = "MS_code_for_neandertal_from_Fabrizio.txt"), collapse = "")

ms <- regmatches(MS_Code,regexec(pattern = "ms (.*?) -t",MS_Code))
ms <- strsplit(ms[[1]][2],split = " ")
ms <- cbind(ms[[1]][1],ms[[1]][2])
ms <- as.numeric(ms)

Theta <-  regmatches(MS_Code,regexec(pattern = "-t (.*?) -r",MS_Code))
Theta <- Theta[[1]][2]
Theta <- as.numeric(Theta)

Recomb_Rate <- regmatches(MS_Code,regexec(pattern = "-r (.*?) -en",MS_Code))
Recomb_Rate <- strsplit(Recomb_Rate[[1]][2],split = " ")
Recomb_Rate <- cbind(Recomb_Rate[[1]][1],Recomb_Rate[[1]][2])
nSites <- as.numeric(Recomb_Rate[2])

N0=Theta/(4*Mut_rate_Vindjia*nSites)

Recomb_Rate <- as.numeric(Recomb_Rate[1])/(4*N0*nSites)

Get_Demographic_Event_Population_Change_F=function(MS_Code,Theta,N0){
  Demographic_Event <- str_match_all(string = MS_Code,pattern = "-en (.*?) ")
  Time <- as.numeric(Demographic_Event[[1]][,2])
  Demographic_Event <- as.character(Demographic_Event[[1]][,2])
  Demographic_Event <- as.data.frame(Demographic_Event)
  Subpop <- c()
  for (i in Demographic_Event[,1]){
    i=as.character(i)
    x <- str_match_all(string = MS_Code,paste(i," (.*?) ",sep = ""))
    Subpop = c(Subpop, as.numeric(x[[1]][2]))
  }
  Demographic_Event <- cbind(Demographic_Event,Subpop)
  SubPopSize <- c()
  for (i in Demographic_Event[,1]){
    j=Demographic_Event$Subpop[Demographic_Event$Demographic_Event == i]
    x <- str_match_all(string = MS_Code,paste(i,j,"(.*?) ",sep = " "))
    SubPopSize = c(SubPopSize, as.numeric(x[[1]][2]))
  }
  
  Demographic_Event_Population_Change <- as.data.frame(cbind(Time,Subpop,SubPopSize))
  
  Demographic_Event_Population_Change$Time=as.numeric(Demographic_Event_Population_Change$Time)*(4*N0)
  Demographic_Event_Population_Change$SubPopSize=Demographic_Event_Population_Change$SubPopSize*N0
  return(Demographic_Event_Population_Change)
}


Demographic_Events_Vindjia=Get_Demographic_Event_Population_Change_F(MS_Code,Theta,N0)
Vindjia_Demography_unbiased_from_Mut_Rate <- as.data.frame(cbind(Demographic_Events_Vindjia[,1]*Mut_rate_Vindjia,Demographic_Events_Vindjia[,3]*Mut_rate_Vindjia))
names(Vindjia_Demography_unbiased_from_Mut_Rate) <- c("time","PopSize")

write.table(Demographic_Events_Vindjia,file = "~/Desktop/Neandertal_Human_Introgression_Project/Introgression_SIM_and_Analysis/Snakemake_SIM/Inferred_Pop_History/Vindjia_PSMC_inferred_Demography.txt")

write.table(Vindjia_Demography_unbiased_from_Mut_Rate,file = "~/Desktop/Neandertal_Human_Introgression_Project/Introgression_SIM_and_Analysis/Snakemake_SIM/Inferred_Pop_History/Vindjia_PSMC_inferred_Demography_unbiased_from_Mut_Rate.txt")

Vindjia_Demography_unbiased_from_Mut_Rate[,1]/2e-08




