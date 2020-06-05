#Snakemake called programm to read interpolate simulated SNPs on CEU specific recomb map in HapMap format and create .snp file
args <- commandArgs(TRUE)
Recomb_Map_Path <- as.character(args[1])
SNP_Map_Path <- as.character(args[2]) 		
Output <- as.character(args[3])  		

Recomb_Map <- read.table(Recomb_Map_Path,header=T)
SNP_Map <- read.table(SNP_Map_Path,header=F)

SNP_Map_corrected <- c()
for(chr in 1:22){
  Recomb_Map_chr <- Recomb_Map[Recomb_Map$Chr==paste("chr",chr,sep=""),]
  SNP_Map_chr <- SNP_Map[SNP_Map$V2==chr,]
  SNP_Map_chr$V3 <- approx(x=Recomb_Map_chr$Position.bp., y=Recomb_Map_chr$Map.cM., xout=SNP_Map_chr$V4,method = "constant")[[2]]/100
  SNP_Map_corrected <- rbind(SNP_Map_corrected,SNP_Map_chr)
  
}

write.table(x=SNP_Map_corrected,file = Output,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
