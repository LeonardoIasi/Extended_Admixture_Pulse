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

#' Fitting the simple to segment length data
#'
#' This function takes the length of segments and estimate the time since the gene flow 
#' under the assumption of a one generation pulse.
#'
#' @param Segments file at least one column containing the length of each unique segment in cM named [length_cM]
#' @param lower_trunc optional parameter to define lower cutoff for the segment length given in cM
#' @param upper_trunc optional parameter to define upper cutoff for the segment length given in cM
#' @param tm_lower gives the lowest values the optimization considers for the mean admixture time
#' @param tm_upper gives the highest values the optimization considers for the mean admixture time
#' @param k_lower gives the lowest values the optimization considers for the shape parameter k
#' @param k_upper gives the highest values the optimization considers for the shape parameter k
#' @examples
#' Segments <- read.table("Example_seg.txt",stringsAsFactors = F,header=T,fill = T)
#' fit_simple_pulse(Segments=Segments,truncation=T,lower_trunc=0.05,lower_trunc=1.2,tm_lower=100,tm_upper=5000)


fit_simple_pulse=function(Segments,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper){
  
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

#' Fitting the extended pulse to segment length data
#'
#' This function takes the length of segments and estimate the time since the gene flow 
#' with gene flow potentially between 1 generation and a constant gene flow
#'
#' @param Segments file at least one column containing the length of each unique segment in cM named [length_cM]
#' @param lower_trunc optional parameter to define lower cutoff for the segment length given in cM
#' @param upper_trunc optional parameter to define upper cutoff for the segment length given in cM
#' @param tm_lower gives the lowest values the optimization considers for the mean admixture time
#' @param tm_upper gives the highest values the optimization considers for the mean admixture time
#' @param k_lower gives the lowest values the optimization considers for the shape parameter k
#' @param k_upper gives the highest values the optimization considers for the shape parameter k
#' Segments <- read.table("Example_seg.txt",stringsAsFactors = F,header=T,fill = T)
#' fit_extended_pulse(Segments=Segments,truncation=T,lower_trunc=0.05,lower_trunc=1.2,tm_lower=100,tm_upper=5000,k_lower=2,k_upper=1e10)

fit_extended_pulse=function(Segments,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper,k_lower,k_upper){
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
    
    res = optim(par = c(500,(1/20)), f, method="L-BFGS-B",lower = c(tm_lower,1/k_upper),upper = c(tm_lower,1/k_lower))
    tm = res$par[1]
    k = 1/res$par[2]
    return(list(res,data.frame(k=k,tm=tm, ll_lomax = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc)))
  }, silent=F)
  return(list(NA,data.frame(k=NA,tm=NA, ll_lomax=NA,lower_trunc=lower_trunc,upper_trunc=upper_trunc)))
}
