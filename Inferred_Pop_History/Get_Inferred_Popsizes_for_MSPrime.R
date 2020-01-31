### MSMC inferred Popsizes of YRI and CEU converted to MSprime Code ###
#install.packages("xlsx")
library("xlsx")
Generation <- 30
Mutation_Rate_Humans <- 1.25e-8

All_inferred_Pops <- read.xlsx(file = "Inferred_Pop_History/MSMC_Inferred_Pop_Sizes.xlsx",sheetIndex = 1,colClasses = c("character","character","integer","character","numeric","numeric","numeric","numeric","logical"))

Yorubans_In <- All_inferred_Pops[41:80,]
CEU_In <- All_inferred_Pops[161:200,]

CEU_time <- CEU_In[,6]/30
CEU_Out <- as.data.frame(cbind(CEU_time,CEU_In[,8]))
names(CEU_Out) <- c("time","PopSize")
write.table(file = "Inferred_Pop_History/CEU_MSMC_inferred_Demography.txt",x=CEU_Out)
CEU_Unbiased_from_Mut_Rate <- as.data.frame(cbind(CEU_Out[,1]*Mutation_Rate_Humans,CEU_Out[,2]*Mutation_Rate_Humans))
write.table(file = "Inferred_Pop_History/CEU_MSMC_inferred_Demography_unbiased_from_Mut_Rate.txt",x=CEU_Unbiased_from_Mut_Rate)

Yorubans_time <- Yorubans_In[,6]/30
Yorubans_Out <- as.data.frame(cbind(Yorubans_time,Yorubans_In[,8]))
names(Yorubans_Out) <- c("time","PopSize")
write.table(file = "Inferred_Pop_History/Yorubans_MSMC_inferred_Demography.txt",x=Yorubans_Out)
Yorubans_Unbiased_from_Mut_Rate <- as.data.frame(cbind(Yorubans_Out[,1]*Mutation_Rate_Humans,Yorubans_Out[,2]*Mutation_Rate_Humans))
write.table(file = "Inferred_Pop_History/Yorubans_MSMC_inferred_Demography_unbiased_from_Mut_Rate.txt",x=Yorubans_Unbiased_from_Mut_Rate)