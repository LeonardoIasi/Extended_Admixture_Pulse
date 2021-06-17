#' Fitting the simple to segment length data
#'
#' This function takes the length of segments and estimate the time since the gene flow 
#' under the assumption of a one generation pulse.
#'
# param: input = Segments filepath with at least one column containing the length of each unique segment in cM named [length_cM]
# param: truncation = optional parameter (TRUE/FALSE) to define if you want to exclude segments given a ceartain threshold lengtn in cM from the fitting (defined yu upper and lower trunc)
# param: lower_trunc = (if truncation = T) numeric parameter to define lower cutoff for the segment length given in cM
# param: upper_trunc = (if truncation = T) numeric parameter to define upper cutoff for the segment length given in cM
# param: tm_lower = gives the lowest values (numeric) the optimization considers for the mean admixture time
# param: tm_upper = gives the highest values (numeric) the optimization considers for the mean admixture time
# param: k_lower = (only for fitting the extended pulse) gives the lowest values (numeric) the optimization considers for the shape parameter k
# param: k_upper = (only for fitting the extended pulse) gives the highest values (numeric) the optimization considers for the shape parameter k

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

input <- "Example_seg.txt"

Output_prefix <- as.character("Example_seg_output")  		
truncation <-  T
lower_trunc <- 0.05	
upper_trunc <- 1
tm_lower <-  100
tm_upper <-  5000
k_lower <- 2
k_upper <- 1e8

fit_simple_pulse=function(Segments,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper){
  Segments <- read.table(Segments,header = T)
  
  if(truncation==T){
    l=Segments$length_cM[Segments$length_cM>=lower_trunc & Segments$length_cM<=upper_trunc]
  } else {
    l=Segments$length_cM[Segments$length_cM>0]
  }
  
  dtexp <- function(x, rate=1, lower_trunc, upper_trunc){
    dexp(x, rate, log=T) - logspace.sub(pexp(lower_trunc, rate, lower=F, log=T),pexp(upper_trunc, rate, lower=F, log=T))}
  
  dtexp_norm <- function(x, rate=1, lower_trunc, upper_trunc){dexp(x, rate, log=T)} 
  
  try({
    
    if(truncation==T){
      f = function(par) -sum(dtexp(l/100, par[1], lower_trunc=lower_trunc/100,upper_trunc = upper_trunc/100))
    } else {
      f = function(par) -sum(dtexp_norm(l/100, par[1]))
    }
    
    res = optim(c(500), f, method="L-BFGS-B",lower = c(tm_lower),upper = c(tm_upper))
    par =res$par
    t_m = par
    return(list(res,data.frame(t_m=t_m, ll_exp = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc)))
  }, silent=F)
  return(list(NA,data.frame(t_m=NA, ll_exp=NA ,lower_trunc=lower_trunc,upper_trunc=upper_trunc)))
}

fit_extended_pulse=function(Segments,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper,k_lower,k_upper){
  Segments <- read.table(Segments,header = T)
  if(truncation==T){
    l=Segments$length_cM[Segments$length_cM>=lower_trunc & Segments$length_cM<=upper_trunc]
  } else {
    l=Segments$length_cM[Segments$length_cM>0]
  }
  
  dtlomax <- function(x, scale, shape, lower_trunc, upper_trunc){
    VGAM::dlomax(x, shape3.q = shape,scale = scale, log=T) -
      logspace.sub(VGAM::plomax(lower_trunc, scale = scale, shape3.q = shape, lower=F, log=T),VGAM::plomax(upper_trunc, scale = scale, shape3.q = shape, lower=F, log=T))}
  
  dtlomax_norm <- function(x, scale, shape)
    VGAM::dlomax(x, shape3.q = shape,scale = scale, log=T) 
  
  
  try({
    if(truncation==T){
      f = function(par)-sum(dtlomax(l/100, scale=((1/par[2])/par[1]), shape=((1/par[2])+1), lower_trunc=lower_trunc/100,upper_trunc=upper_trunc/100))
    } else {
      f = function(par)-sum(dtlomax_norm(l/100, scale=((1/par[2])/par[1]), shape=((1/par[2])+1) ))
    }
    
    res = optim(par = c(500,(1/20)), f, method="L-BFGS-B",lower = c(tm_lower,1/k_upper),upper = c(tm_upper,1/k_lower))
    tm = res$par[1]
    k = 1/res$par[2]
    return(list(res,data.frame(k=k,tm=tm, ll_lomax = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc)))
  }, silent=F)
  return(list(NA,data.frame(k=NA,tm=NA, ll_lomax=NA,lower_trunc=lower_trunc,upper_trunc=upper_trunc)))
}


# Fitting the simple pulse

SP_Fit <- fit_simple_pulse(input,truncation ,lower_trunc,upper_trunc,tm_lower,tm_upper)

tm_SP_est <- SP_Fit[[2]]$t_m
ll_SP_est <- SP_Fit[[2]]$ll_exp
lower_trunc_tm_SP_est <- SP_Fit[[2]]$lower_trunc
upper_trunc_tm_SP_est <- SP_Fit[[2]]$upper_trunc



Fit_simple_pulse_df <- data.frame(tm = tm_SP_est,log_likelihood = ll_SP_est,lower_trunc=lower_trunc_tm_SP_est,upper_trunc=upper_trunc_tm_SP_est)


write.csv(Fit_simple_pulse_df,file = paste(Output_prefix,"_simple_pulse.txt",sep = ""),quote = F,row.names = F)



# Fit the extended pulse

EP_Fit <- fit_extended_pulse(input,truncation,lower_trunc,upper_trunc,tm_lower,tm_upper,k_lower,k_upper)

tm_EP_est <- EP_Fit[[2]]$tm
k_EP_est <- EP_Fit[[2]]$k
ll_EP_est <- EP_Fit[[2]]$ll_lomax
lower_trunc_tm_EP_est <- EP_Fit[[2]]$lower_trunc
upper_trunc_tm_EP_est <- EP_Fit[[2]]$upper_trunc

Fit_extended_pulse_df <- data.frame(tm = tm_EP_est,k = k_EP_est,log_likelihood = ll_EP_est,lower_trunc=lower_trunc_tm_EP_est,upper_trunc=upper_trunc_tm_EP_est)

write.csv(Fit_extended_pulse_df,file = paste(Output_prefix,"_extended_pulse.txt",sep = ""),quote = F,row.names = F)

