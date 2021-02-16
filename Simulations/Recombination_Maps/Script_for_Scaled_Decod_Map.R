# Scale decode Map
startMap=10e6
lengthMap=150e6
Recomb_Map <- read.table("RecombmapIceland.txt", header = T)

Ice_chr_1 <- Recomb_Map[Recomb_Map$Chr=="chr1",]
Decode_Map_chr_1<- data.frame(Ice_chr_1$Chr,Ice_chr_1$Begin,Ice_chr_1$cMperMb,Ice_chr_1$cM)
colnames(Decode_Map_chr_1) <- c("Chromosome","position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")

Cleaved_map <- Decode_Map_chr_1[Decode_Map_chr_1$position>= startMap,]
Starting_Position=Cleaved_map$position[1]
Starting_distance=Cleaved_map$`Genetic_Map(cM)`[1]
Cleaved_map$position=Cleaved_map$position- Starting_Position
Cleaved_map$`Genetic_Map(cM)`=Cleaved_map$`Genetic_Map(cM)`- Starting_distance
Scaled_Map=Cleaved_map[Cleaved_map$position<= lengthMap,]
Scaled_Map$position[1] <- 1
Last_entrie_Genetic_Map=(Scaled_Map$`COMBINED_rate(cM/Mb)`[nrow(Scaled_Map)]* ((lengthMap-Scaled_Map$position[nrow(Scaled_Map)])/1e06) )  + Scaled_Map$`Genetic_Map(cM)`[nrow(Scaled_Map)]
Scaled_Map <- rbind(Scaled_Map,c("chr1",lengthMap,0,Last_entrie_Genetic_Map))
colnames(Scaled_Map)<- c("Chromosome","position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")
write.table(Scaled_Map,"../Recombination_Maps/Scaled_Decode_Map.txt",quote = F,sep = "\t",row.names = F,col.names = T)
                   