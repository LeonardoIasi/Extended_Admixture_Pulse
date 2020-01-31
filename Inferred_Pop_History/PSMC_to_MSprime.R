Length=150e6 
mu=2e-8 
ipop=2 
N0ref=1000
bs=1000
thetapsmc=0.00930765

PSMC_Inferred <- read.table("Desktop/Neandertal_Human_Introgression_Project/Paper/Inferred_Pop_History/Vindija33.19.corrected_read_in.hist")
N0 <- PSMC_Inferred$V3[1]

ThetaMS <- thetapsmc*N0*Length/100

N0psmc=thetapsmc/(4*100*mu)
N0ms=ThetaMS/(4*mu)/Length
lambdarescaling=N0ms/N0ref
time <- bs/(4*N0ref)+PSMC_Inferred$V2[1:length(PSMC_Inferred$V2)]*lambdarescaling*N0psmc/N0ms/2
PopSize <-PSMC_Inferred$V3[1:length(PSMC_Inferred$V3)]*lambdarescaling/N0

MS_Demography <- data.frame(time,PopSize)
MS_Demography$time=as.numeric(MS_Demography$time)*(4*N0ms)
MS_Demography$PopSize=MS_Demography$PopSize*N0ms
MS_Demography_1 <- data.frame(0,MS_Demography$PopSiz[1])
colnames(MS_Demography_1) <- c('time','PopSize')
MS_Demography <- rbind(MS_Demography_1,MS_Demography)

Vindjia_Demography_unbiased_from_Mut_Rate <- as.data.frame(cbind(MS_Demography[,1]*mu,MS_Demography[,2]*mu))
names(Vindjia_Demography_unbiased_from_Mut_Rate) <- c("time","PopSize")

write.table(MS_Demography,file = "~/Desktop/Neandertal_Human_Introgression_Project/Paper/Inferred_Pop_History/Vindjia_PSMC_inferred_Demography_new.txt")

write.table(Vindjia_Demography_unbiased_from_Mut_Rate,file = "~/Desktop/Neandertal_Human_Introgression_Project/Paper/Inferred_Pop_History/Vindjia_PSMC_inferred_Demography_unbiased_from_Mut_Rate_new.txt")


