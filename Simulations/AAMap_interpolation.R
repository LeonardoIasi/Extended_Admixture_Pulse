#Snakemake called programm to read interpolate simulated SNPs on AAMap and create .snp file
args <- commandArgs(TRUE)
#setwd("/home/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Check_Cox_2019_Paper/")
SNP_Map_Path <- as.character(args[1]) 		# output from rollof
Recombination_Map_path <- as.character(args[2])  		# output of expfit log file
Output <- as.character(args[3])  		# output of expfit table of calculated values

SNP_Map <- read.table(file = SNP_Map_Path,header = FALSE)
Recombination_Map <- read.table(file = Recombination_Map_path,header = TRUE)

#library(foreach)
#library(doMC)
#registerDoMC(8)
#SNP_Map[,3] <- foreach(i=SNP_Map[,4], .combine=rbind) %dopar%{
#  
#  xx= tail(Recombination_Map$Map.cM.[Recombination_Map$Position.bp.<=i],n = 1)
#  
#  xx
#}
#SNP_Map[,3] <- Recombination_Map[vapply(SNP_Map[,4], function(x) which(x <= Recombination_Map[, 2])[1] - 1L, 
#                                        FUN.VALUE = integer(1)),4]/100

SNP_Map[,3] <- approx(x=Recombination_Map$Position.bp., y=Recombination_Map$Map.cM., xout=SNP_Map[,4],method = "constant")[[2]]/100

New_SNP_Map <- as.data.frame(SNP_Map)

write.table(x=New_SNP_Map,file = Output,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)


