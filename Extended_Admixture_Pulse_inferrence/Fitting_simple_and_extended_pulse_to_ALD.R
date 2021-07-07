# Fitting the simple to ALD length data
#
# This function takes the Ancestry Linkage Disequilibrium and estimate the time since the gene flow 
# under the assumption of a one generation pulse.
#
# param: input =  ALD filepath with the first column containing the distance of SNP's in cM and a the second column the weighted LD
# param: lower_trunc =  optional parameter to define a lower cutoff for the ALD length given in cM
# param: upper_trunc = optional parameter to define upper cutoff for the ALD length given in cM
# param: constant = optional parameter defines if a constant c is fitted to model background LD (recommended)
# param: tm_lower = gives the lowest values the optimization considers for the mean admixture time
# param: tm_upper = gives the highest values the optimization considers for the mean admixture time
# param: k_lower = (only for fitting the extended pulse) gives the lowest values the optimization considers for the shape parameter k
# param: k_upper  = (only for fitting the extended pulse) gives the highest values the optimization considers for the shape parameter k

suppressPackageStartupMessages({
  library(VGAM)
  library(tidyverse)
  library("DEoptim")
  library("MASS")
  library(bbmle)
  library(rethinking)
  library(DPQ)
  library(viridis)
})

Get_points <- function(input,lower_trunc,upper_trunc){
  # Read input file
  data <- read.table(input, header = F)
  
  
  # set dist and LD
  dist <- data[,1]
  LD <- data[,2]
  ndist <- length(dist)  ## number of rows in dataset
  lower_trunc=lower_trunc
  upper_trunc=upper_trunc
  
  # check x lower value and y lower value
  data.sub <- data
  if ((lower_trunc > dist[1]) || (upper_trunc < dist[ndist])) {
    data.sub <- subset(data, ((dist <= upper_trunc) & (dist >= lower_trunc)))
  }
  dist <- data.sub[,1]		# updated x values
  LD <- data.sub[,2]		# updated y values
  result_table <- data.frame(dist,LD)
  result_table <- result_table[complete.cases(result_table),]
  return(result_table)
}


Simple_Pulse_ALD_fn <- function(Data,constant,tm_lower,tm_upper){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$LD
  
  if(constant){
    fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])+x[3]
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,tm_lower,-1), upper=c(1, tm_upper,1), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    names(par1_exp) <- c("A", "tm","c")
    
    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)+c), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE)) 
    f_predict = data.frame(dist_M = xx$dist/100, ALD = predict(fit1_exp))
  } else {
    fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,tm_lower), upper=c(1, tm_upper), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    names(par1_exp) <- c("A", "tm")
    
    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))
    f_predict = data.frame(dist_M = xx$dist/100, ALD = predict(fit1_exp))
    
  }
  return(list(fit1_exp,data.frame(tm=coef(fit1_exp)[2], RSS = sum(resid(fit1_exp)^2),lower_trunc=min(xx$dist),upper_trunc=max(xx$dist)),xx,f_predict))
}

Extended_Pulse_ALD_fn <- function(Data,constant,tm_lower,tm_upper,k_lower,k_upper){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$LD
  
  if(constant){
    fm1_lomax <- function(x) x[4] + x[3]* (1/(1 + (  (x[2] / (1/x[1]))*(dist/100) )))^((1/x[1]))
    fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1/k_upper,tm_lower,1e-6,-1), upper=c(1/k_lower,tm_upper,1,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],par1_lomax[2],par1_lomax[3],par1_lomax[4])
    names(par1_lomax) <- c("k","tm","A","c")
    fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((tm/ (1/k)) *(dist/100)) ))^((1/k)), start=par1_lomax, algorithm="port",
                      lower=c(1/k_upper,tm_lower,1e-6,-1), upper=c(1/k_lower,tm_upper,1,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
    f_predict = data.frame(dist_M = xx$dist/100, ALD = predict(fit1_lomax))
    
  } else {
    fm1_lomax <- function(x)  x[3]* (1/(1 + (  (x[2] / (1/x[1]))*(dist/100) )))^((1/x[1]))
    fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1/k_upper,tm_lower,1e-6,-1), upper=c(1/k_lower,tm_upper,1,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],par1_lomax[2],par1_lomax[3])
    names(par1_lomax) <- c("k","tm","A")
    fit1_lomax <- nls(wcorr ~ A*(1/(1  + ((tm/ (1/k)) *(dist/100)) ))^((1/k)), start=par1_lomax, algorithm="port",
                      lower=c(1/k_upper,tm_lower,1e-6,-1), upper=c(1/k_lower,tm_upper,1,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
    f_predict = data.frame(dist_M = xx$dist/100, ALD = predict(fit1_lomax))
    
  }
  return(list(fit1_lomax,data.frame(k=1/coef(fit1_lomax)[1],tm=coef(fit1_lomax)[2],RSS = sum(resid(fit1_lomax)^2),lower_trunc=min(xx$dist),upper_trunc=max(xx$dist)),xx,f_predict))
}


Get_CI_ALD_fn <- function(est,n_data){
  org <- est
  lwr_approx <- est*(1-(1.96/sqrt(length(n_data[,1]))))
  upr_approx <- est*(1+(1.96/sqrt(length(n_data[,1]))))
  param_CI <- c(org,lwr_approx,upr_approx)
  return(param_CI)
}
