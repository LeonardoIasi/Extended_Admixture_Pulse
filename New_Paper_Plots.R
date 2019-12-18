#### New Figure for Paper using recent sampling (50 gen): constant RR, RM no correction, RM diff RM correction, RM correction same RM ####

setwd('Desktop/Neandertal_Human_Introgression_Project/Paper/')
library(ggplot2)
library(reshape2)
library(viridis)
library(ggpubr)
cbPalette_viridis <- viridis(6,option = "D")

Filter_Result_Table <- function(Result.Table){
  Result.Table=read.table(Result.Table,header=F,sep = " ")
  Result.names <- c('A.exp', 's.exp', 'c.exp', 'RSS_Expo','AIC_Expo', 'A.lomax', 's.lomax', 'w.lomax','c.lomax', 'RSS_Lomax','AIC_Lomax', 'F_Test','Scenario.name','GF.Start','GF.End','AS','minDist','GF.Model')
  colnames(Result.Table) <- Result.names
  #Result.Table <- Result.Table[(1/Result.Table$s.exp)>0,]
  #Result.Table <- Result.Table[(1/Result.Table$s.exp)<5000,]
  #Result.Table <- Result.Table[Result.Table$RSS_Lomax<1e-2,]
  return(Result.Table)
}

#### True values
True_params <- function(Result.Table){
  True_params_Calc <- function(True_GF_length,True_mean_GF){
    EX= True_mean_GF
    GF_Len <-True_GF_length
    VarX= ((GF_Len)/4)**2
    b= EX/VarX
    a=EX*b
    #a=a+1
    True_W=1/a
    True_S=b/(1/True_W)
    xx=c(True_W,True_S)
    return(xx)
  }
  True.params <- c()
  for (i in 1:length(Result.Table$F_Test)) {
    xx=True_params_Calc(Result.Table$True_GF_length[i],Result.Table$True_mean_GF[i])
    True.params <- rbind(True.params,xx)
  }
  True.params <- as.data.frame(True.params)
  return(True.params)
}

Result.Table.fn <- function(Result.Table.path,Sampling.time.from.GF.End,name){
  Result.Table <- Filter_Result_Table(Result.Table.path)
  Result.Table$Name <- name
  Result.Table$Sample_Time <- Sampling.time.from.GF.End
  Result.Table$True_GF_length <- Result.Table$GF.End-Result.Table$GF.Start
  Result.Table$True_mean_GF <- ((Result.Table$GF.End+Result.Table$GF.Start)/2)-(Result.Table$GF.Start-Result.Table$Sample_Time)
  Result.Table$mean_GF_exp <- 1/Result.Table$s.exp
  Result.Table$mean_GF_lomax_s <- 1/Result.Table$s.lomax
  Result.Table.comparison <- True_params(Result.Table)
  Result.Table$True_W <- Result.Table.comparison$V1
  Result.Table$True_S <- Result.Table.comparison$V2
  Result.Table.comparison$length_GF <- (sqrt(Result.Table$mean_GF_lomax_s/((1/Result.Table$w.lomax)*Result.Table$s.lomax)))*4
  Result.Table$length_GF <- Result.Table.comparison$length_GF
  
  return(Result.Table)
}

Plot_sampling_50_gen_t_GF <- function(ggdata_t.GF,cbPalette_viridis){
  ggplot(ggdata_t.GF,aes(x=variable,y=as.numeric(as.character(value)),colour=factor(Name)))+
    geom_boxplot(show.legend = T)+
    facet_grid(~as.factor(True_GF_length),switch = "x")+
    geom_hline( aes(yintercept = True_mean_GF ))+
    ggtitle("All sampled 50 gen from GF end") +
    labs(x = "Gene Flow Length")+
    labs(y = "Estimated Admixture Time")+
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    coord_cartesian(ylim = c(0,1000), expand = 0)+
    scale_color_manual("Simulations",
                       values = c(cbPalette_viridis[1],cbPalette_viridis[2],cbPalette_viridis[3],cbPalette_viridis[4]),
                       labels = c("Constant Recombination","HapMap AAMap corrected","HapMap HapMap corrected","HapMap not corrected"))
  
}

Plot_sampling_50_gen_GF_Length <- function(ggdata_l.GF,cbPalette_viridis){
  ggplot(ggdata_l.GF,aes(x=variable,y=as.numeric(as.character(value)),colour=factor(Name)))+
    geom_boxplot(show.legend = F)+
    facet_grid(~as.factor(True_GF_length),switch = "x")+
    geom_hline( aes(yintercept = True_GF_length ))+
    ggtitle(paste("All sampled 50 gen from GF end",sep = " ")) +
    labs(x = "Gene Flow Length")+
    labs(y = "Estimated Gene Flow Length")+
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    coord_cartesian(ylim = c(0,4000), expand = 0)+
    scale_color_manual("Simulations",
                       values = c(cbPalette_viridis[1],cbPalette_viridis[2],cbPalette_viridis[3],cbPalette_viridis[4]),
                       labels = c("Constant Recombination","HapMap AAMap corrected","HapMap HapMap corrected","HapMap not corrected"))
}

Plot_sampling_varying_gen_t_GF <- function(ggdata_t.GF_varying,cbPalette_viridis){
  ggplot(ggdata_t.GF_varying,aes(x=variable,y=as.numeric(as.character(value)),colour=factor(Name)))+
    geom_boxplot()+
    facet_grid(~as.factor(Sample_Time),switch = "x")+
    geom_hline( aes(yintercept = True_mean_GF ))+
    ggtitle("All Gene Flow length = 800 gen") +
    labs(x = "Generations Sample After Gene Flow End")+
    labs(y = "Estimated Admixture Time")+
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    coord_cartesian(ylim = c(0,2000), expand = 0)+
    scale_color_manual("Simulations",
                       values = c(cbPalette_viridis[1],cbPalette_viridis[2],cbPalette_viridis[3],cbPalette_viridis[4]),
                       labels = c("Constant Recombination","HapMap AAMap corrected","HapMap HapMap corrected","HapMap not corrected"))
  
}

Plot_sampling_varying_gen_GF_Length <- function(ggdata_l.GF_varying,cbPalette_viridis){
  ggplot(ggdata_l.GF_varying,aes(x=variable,y=as.numeric(as.character(value)),colour=factor(Name)))+
    #geom_point(aes(size = Freq),show.legend = F)+
    geom_boxplot(show.legend = F)+
    facet_grid(~as.factor(Sample_Time),switch = "x")+
    geom_hline( aes(yintercept = True_GF_length ))+
    ggtitle("All Gene Flow length = 800 gen") +
    labs(x = "Generations Sample After Gene Flow End")+
    labs(y = "Estimated Gene Flow  Length")+
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    coord_cartesian(ylim = c(0,4000), expand = 0)+
    scale_color_manual("Simulations",
                       values = c(cbPalette_viridis[1],cbPalette_viridis[2],cbPalette_viridis[3],cbPalette_viridis[4]),
                       labels = c("Constant Recombination","HapMap AAMap corrected","HapMap HapMap corrected","HapMap not corrected"))
}

Result.Table.path_Recent_50 <-'Close_to_GF_End_Recent_GF/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Close_to_GF_End_Recent_GF-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-No_correction.txt'
Result.Table.path_Recent_50_Recomb_Map_no_correction <-'Close_to_GF_End_Recent_GF_Recomn_Map_Hap_Map/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Close_to_GF_End_Recent_GF_Recomn_Map_Hap_Map-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-No_correction.txt'
Result.Table.path_Recent_50_Recomb_Map_AAMap_correction <-'Close_to_GF_End_Recent_GF_Recomn_Map_Hap_Map/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Close_to_GF_End_Recent_GF_Recomn_Map_Hap_Map-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-AAMap_correction.txt'
Result.Table.path_Recent_50_Recomb_Map_HapMap_Recomb_HapMap_correction <-'Close_to_GF_End_Recent_GF_Recomn_Map_Hap_Map/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Close_to_GF_End_Recent_GF_Recomn_Map_Hap_Map-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-HapMap_correction.txt'

Result.Table.path_Recent_Varying <-'Variyng_Time_of_Recent_sampling_GF_l_fixed/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Variyng_Time_of_Recent_sampling_GF_l_fixed-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-No_correction.txt'
Result.Table.path_Recent_Varying_Recomb_Map_no_correction <-'Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-No_correction.txt'
Result.Table.path_Recent_Varying_Recomb_Map_AAMap_correction <-'Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-AAMap_correction.txt'
Result.Table.path_Recent_Varying_Recomb_Map_HapMap_Recomb_HapMap_correction <-'Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map/Result_both_Fit/Result_file_SIM_Raw_ALDER-Fit-Variyng_Time_of_Recent_sampling_GF_l_fixed_Recomn_Map_Hap_Map-GF_Model_IV-min_dist_Fit-0.05-ascertainment-0-HapMap_correction.txt'

Plot.data_50<-rbind(
  Recent_50_constant<- Result.Table.fn(Result.Table.path_Recent_50,50,'Recent_50_constant'),
  Recent_50_HapMap_no_correction <- Result.Table.fn(Result.Table.path_Recent_50_Recomb_Map_no_correction,50,'Recent_50_HapMap_no_correction'),
  Recent_50_HapMap_AAMap_correction <- Result.Table.fn(Result.Table.path_Recent_50_Recomb_Map_AAMap_correction,50,'Recent_50_HapMap_AAMap_correction'),
  Recent_50_HapMap_HapMap_correction <- Result.Table.fn(Result.Table.path_Recent_50_Recomb_Map_HapMap_Recomb_HapMap_correction,50,'Recent_50_HapMap_HapMap_correction')
)

ggdata_t.GF_50 <- melt(Plot.data_50,measure.vars =  c('mean_GF_exp'),id.vars = c('True_mean_GF','True_GF_length','Name'))
ggdata_l.GF_50 <- melt(Plot.data_50,measure.vars =  c('length_GF'),id.vars = c('True_mean_GF','True_GF_length','Name'))

Plot.data_varying<-rbind(
  Recent_varying_constant<- Result.Table.fn(Result.Table.path_Recent_Varying,c(rep(50,100),rep(100,100),rep(200,100),rep(400,100),rep(800,100),rep(1000,100)),'Recent_50_constant'),
  Recent_varying_HapMap_no_correction <- Result.Table.fn(Result.Table.path_Recent_Varying_Recomb_Map_no_correction,c(rep(50,100),rep(100,100),rep(200,100),rep(400,100),rep(800,100),rep(1000,100)),'Recent_50_HapMap_no_correction'),
  Recent_varying_HapMap_AAMap_correction <- Result.Table.fn(Result.Table.path_Recent_Varying_Recomb_Map_AAMap_correction,c(rep(50,100),rep(100,100),rep(200,100),rep(400,100),rep(800,100),rep(1000,100)),'Recent_50_HapMap_AAMap_correction'),
  Recent_varying_HapMap_HapMap_correction <- Result.Table.fn(Result.Table.path_Recent_Varying_Recomb_Map_HapMap_Recomb_HapMap_correction,c(rep(50,100),rep(100,100),rep(200,100),rep(400,100),rep(800,100),rep(1000,100)),'Recent_50_HapMap_HapMap_correction')
)

ggdata_t.GF_varying <- melt(Plot.data_varying,measure.vars =  c('mean_GF_exp'),id.vars = c('True_mean_GF','True_GF_length','Name','Sample_Time'))
ggdata_l.GF_varying <- melt(Plot.data_varying,measure.vars =  c('length_GF'),id.vars = c('True_mean_GF','True_GF_length','Name','Sample_Time'))

New.P1 <- Plot_sampling_50_gen_t_GF(ggdata_t.GF_50,cbPalette_viridis)
New.P2 <- Plot_sampling_50_gen_GF_Length(ggdata_l.GF_50,cbPalette_viridis)
New.P3 <- Plot_sampling_varying_gen_t_GF(ggdata_t.GF_varying,cbPalette_viridis)
New.P4 <- Plot_sampling_varying_gen_GF_Length(ggdata_l.GF_varying,cbPalette_viridis)


New.Fig <- ggarrange(New.P1,New.P2,New.P3,New.P4,
                   labels = c("A","B","C","D"),
                   ncol = 2, nrow = 2,common.legend = T,legend = 'bottom')


pdf('~/Desktop/Diff_Sampling_ATE/All_plot_together.pdf',width = 12)
annotate_figure(New.Fig,
                top = text_grob("All Recomb Sim", color = "black", face = "bold", size = 10))
dev.off()