# Fitting the simple and extended pulse to segment length data
# data (Fragments) must be at least one column containing the length of each unique segment in cM

Segments =  read.table("Your_segment_table.txt",header = T) 
lower_trunc = as.numerical(0.05) # the smallest fragment length that can be reliably detected
upper_trunc = as.numeric(1) # the upper length that can be reliably detected

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


fit_simple_pulse=function(Segments,lower_trunc,upper_trunc){
  l=Fragments[Fragments>=lower_trunc & Fragments<=upper_trunc]
  n = rep(1, length(l))
  
  
  dtexp <- function(x, rate=1, lower_trunc, upper_trunc,m,Ne){
    dexp(x, rate, log=T) - logspace.sub(pexp(lower_trunc, rate, lower=F, log=T),pexp(upper_trunc, rate, lower=F, log=T))
  }

  
  try({
    f = function(par) -sum(n * dtexp(l, exp(par[1]), lower_trunc=lower_trunc,upper_trunc = upper_trunc))
    res = optim(c(100), f, method="L-BFGS-B")
    par =exp(res$par)
    t_m = (par*100)
    return(data.frame(rate=par, ll_exp = -res$value, t_m= t_m))
  }, silent=F)
  return(data.frame(rate=NA, ll_exp=NA, t_m= NA ))
}


fit_extended_pulse=function(Segments,lower_trunc,upper_trunc){
  l=Fragments[Fragments>=lower_trunc & Fragments<=upper_trunc]
  n = rep(1, length(l))
  
  
  #' truncated lomax density
  dtlomax <- function(x, scale, shape, lower_trunc, upper_trunc){
    VGAM::dlomax(x, scale=scale, shape=shape, log=T) - 
      logspace.sub(VGAM::plomax(lower_trunc, scale=scale, shape=shape, lower=F, log=T),VGAM::plomax(upper_trunc, scale=scale, shape=shape, lower=F, log=T))
  }
  
  
  try({
    f = function(par)-sum(n * dtlomax(l, scale=exp(par[1]), shape=exp(par[2]), lower_trunc=lower_trunc,upper_trunc=upper_trunc))
    res = optim(c(100,10), f, method="L-BFGS-B")
    pars = exp(res$par)
    return(data.frame(shape=pars[2],scale=pars[1], ll_lomax = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc))
  }, silent=F)
  return(data.frame(shape=NA,scale=NA, ll_lomax=NA))
}


fit_simple_pulse(data_x$length_cM,lower_trunc,upper_trunc)

fit_extended_pulse=function(data_x$length_cM,lower_trunc,upper_trunc)