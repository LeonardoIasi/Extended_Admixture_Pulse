Path_to_Fragment_File = as.character(snakemake@input[[1]])
replication=as.integer(snakemake@wildcards[['number']])
header = as.logic(snakemake@params[['header']])
recomb_rate = as.numeric(snakemake@params[['recomb_rate']])
Path_to_Recomb_Map= as.character(snakemake@params[['Recombination_Map']])
upper_trunc = as.numeric(snakemake@params[['upper_trunc']])
n_Haploid= as.integer(snakemake@params[['n_Haploids']])
Snake_output=as.character(snakemake@output[[1]])


### Functions
Get_Filtered_Simulated_Fragments_fn <- function(Path_to_Fragment_File,header,recomb_rate,n_Haploid){
  Segments <- read.table(Path_to_Fragment_File,stringsAsFactors = F,header=header)
  if(header==F){
    colnames(Segments) = c('Ind','chrom', 'start','end')
  }
  
  Segments$length_bp <- as.numeric(Segments$end) - as.numeric(Segments$start)
  m=sum(Segments$length_bp)/(length(unique(Segments$chrom))*n_Haploid*150e6)
  Segments$length_cM <- Segments$length_bp*(recomb_rate)*100
  
  Segments_filtered <- c()
  for(chr in 1:22){
    Segments_filtered_chr <- Segments[Segments$chrom==paste("chr",chr,sep = "") & !duplicated(Segments$start) & !duplicated(Segments$end),]
    Segments_filtered <- rbind(Segments_filtered,Segments_filtered_chr)
    rm(Segments_filtered_chr)
  }
  return(list(Segments_filtered=Segments_filtered,m=m))
  
}

Assign_Genetic_distance_Simulation_interpolate_fn <- function(Path_to_Fragment_File,header,Path_to_Recomb_Map,n_Haploid){
  Fragments <- read.table(Path_to_Fragment_File,stringsAsFactors = F,header=header,fill = T)
  if(header==F){
    colnames(Fragments) = c('Ind','chrom', 'start','end')
  }
  Fragment_bp=Fragments$end - Fragments$start
  m=sum(Fragment_bp)/(length(unique(Fragments$chrom))*n_Haploid*150e6)
  Recomb_Map <- read.table(Path_to_Recomb_Map, header = T)
  Fragments_cM <- c()
  for(chr in 1:length(unique(Fragments$chrom))){
    Recomb_Map_chr <- Recomb_Map
    Fragments_chr <- Fragments[Fragments$chrom==paste("chr",chr,sep=""),]
    #Start_Arch_Freg <- approx(x=Recomb_Map_chr[,2]  , y=Recomb_Map_chr[,4], xout=Fragments_chr$start,method = "linear")
    #Stop_Arch_Freg <- approx(x=Recomb_Map_chr[,2], y=Recomb_Map_chr[,4], xout=Fragments_chr$end,method = "linear")
    Start_Arch_Freg <- approx(x=Recomb_Map_chr[,2]  , y=Recomb_Map_chr[,4], xout=Fragments_chr$start,method = "constant")
    Stop_Arch_Freg <- approx(x=Recomb_Map_chr[,2], y=Recomb_Map_chr[,4], xout=Fragments_chr$end,method = "constant")
    Chr_x <- data.frame(chrom=paste("chr",chr,sep = ""),Start_cM=Start_Arch_Freg$y,End_cM=Stop_Arch_Freg$y,length_bp=Fragments_chr$end - Fragments_chr$start)
    Fragments_cM <- rbind(Fragments_cM,Chr_x)
  }
  Fragments_cM_clean <- Fragments_cM[complete.cases(Fragments_cM$Start_cM),]
  Fragments_cM_clean <- Fragments_cM_clean[complete.cases(Fragments_cM_clean$End_cM),]
  Fragments_cM_clean$length_cM <- Fragments_cM_clean$End_cM - Fragments_cM_clean$Start_cM
  Fragments_cM_clean_filtered <- c()
  for(chr in 1:22){
    Fragments_cM_clean_filtered_chr <- Fragments_cM_clean[Fragments_cM_clean$chrom==paste("chr",chr,sep = "") & !duplicated(Fragments_cM_clean$Start_cM) & !duplicated(Fragments_cM_clean$End_cM),]
    Fragments_cM_clean_filtered <- rbind(Fragments_cM_clean_filtered,Fragments_cM_clean_filtered_chr)
    rm(Fragments_cM_clean_filtered_chr)
  }
  return(list(Fragments_cM_clean_filtered=Fragments_cM_clean_filtered,m=m))
  
}

#truncated expo density
dtexp <- function(x, rate=1, lower_trunc, upper_trunc){
  dexp(x, rate, log=T) - logspace.sub(pexp(lower_trunc, rate, lower=F, log=T),pexp(upper_trunc, rate, lower=F, log=T))
}

#truncated lomax density
dtlomax <- function(x, scale, shape, lower_trunc, upper_trunc){
  VGAM::dlomax(x, scale=scale, shape=shape, log=T) - 
    logspace.sub(VGAM::plomax(lower_trunc, scale=scale, shape=shape, lower=F, log=T),VGAM::plomax(upper_trunc, scale=scale, shape=shape, lower=F, log=T))
}

fit_exp_optim_rep=function(Fragments,lower_trunc,upper_trunc){
  l=Fragments[Fragments>=lower_trunc & Fragments<=upper_trunc]
  n = rep(1, length(l))
  
  try({
    f = function(par) -sum(n * dtexp(l, exp(par[1]), lower_trunc=lower_trunc,upper_trunc = upper_trunc))
    #f2 = function(par) -sum(n * dtexp(l, (par[1])))
    res = optim(c(0), f, method="L-BFGS-B")
    par =exp(res$par)
    return(data.frame(rate=par, ll_exp = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc))
  }, silent=F)
  return(data.frame(rate=NA, ll_exp=NA,lower_trunc=lower_trunc,upper_trunc=upper_trunc))
}



fit_lomax_optim_rep=function(Fragments,lower_trunc,upper_trunc){
  l=Fragments[Fragments>=lower_trunc & Fragments<=upper_trunc]
  n = rep(1, length(l))
  
  try({
    f = function(par)-sum(n * dtlomax(l, scale=par[1], shape=par[2],lower_trunc=lower_trunc,upper_trunc=upper_trunc))
    res = optim(c(1, 1), f, method="L-BFGS-B", lower=c(1e-5, 1))
    pars = res$par
    return(data.frame(scale=pars[1], shape=pars[2],ll_lomax=-res$value))
  }, silent=F)
  return(data.frame(scale=NA, shape=NA,ll_lomax=NA))
}

fit_lomax_optim_Kozubowski_rep=function(Fragments,lower_trunc,upper_trunc){
  l=Fragments[Fragments>=lower_trunc & Fragments<=upper_trunc]
  n = rep(1, length(l))
  
  dlomax_Kozubowski <- function(x, scale=scale, shape=shape) {
    s=scale/shape
    w=1/shape
    
    L=1/s * (1/(1+(w*x)/s))^(1/w+1)
    return(L)
    
  }
  
  dlomax_Kozubowski_log <- function(x, scale=scale, shape=shape) {
    s=scale/shape
    w=1/shape
    
    L=(1/w+1) * log(1/(1+ (x*w)/s)) + log(1/s)
    return(L)
    
  }
  
  
  
  plomax_Kozubowski <- function(trunc, scale=scale, shape=shape) {
    s=scale/shape
    w=1/shape
    
    L=(1/(1+(w*trunc)/s))^(1/w)
    return(L)
    
  }
  
  plomax_Kozubowski_log <- function(trunc, scale=scale, shape=shape) {
    s=scale/shape
    w=1/shape
    L=(1/w) * log(1/(1+ trunc*w/s))
    return(L)
    
  }
  #' truncated lomax density
  dtlomax_Kozubowski <- function(x, scale, shape, lower_trunc, upper_trunc){
    dlomax_Kozubowski_log(x, scale=scale, shape=shape) -
      logspace.sub(plomax_Kozubowski_log(lower_trunc, scale=scale, shape=shape),plomax_Kozubowski_log(upper_trunc, scale=scale, shape=shape))
  }
  
  
  try({
    f = function(par)-sum(n * dtlomax_Kozubowski(l, scale=par[1], shape=par[2],lower_trunc=lower_trunc,upper_trunc=upper_trunc))
    res = optim(c(1, 1), f, method="L-BFGS-B", lower=c(1e-8, 1))
    pars = res$par
    return(data.frame(scale_Kozubowski=pars[1], shape_Kozubowski=pars[2],ll_lomax_Kozubowski=-res$value))
  }, silent=F)
  return(data.frame(scale_Kozubowski=NA, shape_Kozubowski=NA,ll_lomax_Kozubowski=NA))
}

Get_Expo_est_fn <- function(Fragments,lower_truncs,upper_trunc,m){
  Expo_Res <- c()
  for(i in lower_truncs){
    expo_fit_res=fit_exp_optim_rep(Fragments$length_cM,i,upper_trunc)
    Expo_Res <- rbind(Expo_Res,expo_fit_res)
  }
  
  Expo_Res$lower_trunc <- lower_truncs
  Expo_Res$upper_trunc <- upper_trunc
  return(Expo_Res)
  
}


Get_Lomax_est_fn <- function(Fragments,Expo_Res,Expo_Res_col_lower_trunc,Expo_Res_col_upper_trunc){
  Res_lomax <- c()
  for(trunc in 1:length(Expo_Res[,1])){
    Res_Segments_optim_lomax=fit_lomax_optim_rep(Fragments$length_cM,Expo_Res[trunc,Expo_Res_col_lower_trunc],Expo_Res[trunc,Expo_Res_col_upper_trunc])
    Res_lomax <- rbind(Res_lomax,Res_Segments_optim_lomax)
  }
  
  return(Res_lomax)
}

Get_Lomax_Kozubowski_est_fn <- function(Fragments,Expo_Res,Expo_Res_col_lower_trunc,Expo_Res_col_upper_trunc){
  Res_lomax <- c()
  for(trunc in 1:length(Expo_Res[,1])){
    Res_Segments_optim_lomax=fit_lomax_optim_Kozubowski_rep(Fragments$length_cM,Expo_Res[trunc,Expo_Res_col_lower_trunc],Expo_Res[trunc,Expo_Res_col_upper_trunc])
    Res_lomax <- rbind(Res_lomax,Res_Segments_optim_lomax)
  }
  
  return(Res_lomax)
}



Results_fit_Simulation <- function(Path_to_Fragment_File,header,Path_to_Recomb_Map,recomb_rate,upper_trunc,n_Haploid,replication){
  if(Path_to_Recomb_Map == "no_Map"){
    Segments_filtered <- Get_Filtered_Simulated_Fragments_fn(Path_to_Fragment_File,header = header,recomb_rate = recomb_rate,n_Haploid=n_Haploid)
    m <- Segments_filtered$m
    Segments_filtered <- Segments_filtered$Segments_filtered
  } else{
    Segments_filtered <- Assign_Genetic_distance_Simulation_interpolate_fn(Path_to_Fragment_File,Path_to_Recomb_Map,header = header,n_Haploid=n_Haploid)
    m <- Segments_filtered$m
    Segments_filtered <- Segments_filtered$Fragments_cM_clean_filtered
  }
  
  Expo_res <- Get_Expo_est_fn(Fragments = Segments_filtered,lower_truncs = c(0.001,0.002,0.005,0.01,0.05,0.1,0.15,0.2,0.5),upper_trunc = upper_trunc,m=m)
  
  Lomax_res <- Get_Lomax_est_fn(Segments_filtered,Expo_res,3,4)
  
  Lomax_Kozubowski_res <- Get_Lomax_Kozubowski_est_fn(Segments_filtered,Expo_res,3,4)
  
  Output=cbind(Expo_res,Lomax_res,Lomax_Kozubowski_res)
  Output$replication_number <- replication
  return(Output)

  
}

Output <- Results_fit_Simulation(Path_to_Fragment_File,header,Path_to_Recomb_Map,recomb_rate,upper_trunc,n_Haploid,replication )

#print output
outlog <- paste(Snake_output)
write.csv(Output,file = outlog,quote = F,row.names = F)


