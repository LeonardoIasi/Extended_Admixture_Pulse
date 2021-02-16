#Snakemake called programm to read interpolate simulated SNPs on AAMap and create .snp file
args <- commandArgs(TRUE)
#setwd("/home/leonardo_iasi/Desktop/Neandertal_Human_Introgression_Project/Check_Cox_2019_Paper/")
SNP_Map_Path <- as.character(args[1]) 		# output from rollof
Recombination_Map_path <- as.character(args[2])  		# output of expfit log file
Output <- as.character(args[3])  		# output of expfit table of calculated values

SNP_Map <- read.table(file = SNP_Map_Path,header = FALSE)
Recombination_Map <- read.table(file = Recombination_Map_path,header = TRUE)

SNP_Map[,3] <- approx(x=Recombination_Map$position, y=Recombination_Map$Genetic_Map.cM., xout=SNP_Map[,4],method = "constant")[[2]]/100

New_SNP_Map <- as.data.frame(SNP_Map)

write.table(x=New_SNP_Map,file = Output,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
