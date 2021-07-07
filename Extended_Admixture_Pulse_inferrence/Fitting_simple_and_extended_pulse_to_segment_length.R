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

Simple_Pulse_Segments_fn=function(input,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper){
  Segments <- read.table(input,header = T)
  
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
    tm = par
    if(truncation==T){
      l_p <- hist(l,breaks = 500,plot = F)
      f_predict = dtexp(l_p$mids/100, tm, lower_trunc=lower_trunc/100,upper_trunc = upper_trunc/100)
      f_predict = data.frame(seg_len_M = l_p$mids, log_dens = f_predict)
    } else {
      l_p <-  hist(l,breaks = 500,plot = F)
      f_predict = dtexp_norm(l_p$mids/100, tm)
      f_predict = data.frame(seg_len_M = l_p$mids, log_dens = f_predict)
    }
    
    return(list(res,data.frame(tm=tm, ll_sp = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc),l,f_predict))
  }, silent=F)
  return(list(NA,data.frame(tm=NA, ll_sp=NA ,lower_trunc=lower_trunc,upper_trunc=upper_trunc),l,NA))
}



Extended_Pulse_Segments_fn=function(input,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper,k_lower,k_upper){
  Segments <- read.table(input,header = T)
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
    if(truncation==T){
      l_p <- hist(l,breaks = 500,plot = F)
      f_predict = dtlomax(l_p$mids/100, scale=(k/tm), shape=(k+1), lower_trunc=lower_trunc/100,upper_trunc=upper_trunc/100)
      f_predict = data.frame(seg_len_M = l_p$mids, log_dens = f_predict)
    } else {
      l_p <-  hist(l,breaks = 500,plot = F)
      f_predict = dtlomax_norm(l/100, scale=(k/tm), shape=(k+1) )
      f_predict = data.frame(seg_len_M = l_p$mids, log_dens = f_predict)
    }
    return(list(res,data.frame(k=k,tm=tm, ll_ep = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc),l,f_predict))
  }, silent=F)
  return(list(NA,data.frame(k=NA,tm=NA, ll_ep=NA,lower_trunc=lower_trunc,upper_trunc=upper_trunc),l,NA))
}

Get_CI_fn <- function(est,n_data){
  org <- est
  lwr_approx <- est*(1-(1.96/sqrt(length(n_data))))
  upr_approx <- est*(1+(1.96/sqrt(length(n_data))))
  param_CI <- c(org,lwr_approx,upr_approx)
  return(param_CI)
}