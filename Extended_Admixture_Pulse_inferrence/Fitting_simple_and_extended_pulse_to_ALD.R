# Fitting the simple to ALD length data
#
# This function takes the Ancestry Linkage Disequilibrium and estimate the time since the gene flow 
# under the assumption of a one generation pulse.
#
# param: input =  ALD filepath with the first column containing the distance of SNP's in cM and a the second column the weighted LD
# param: lval =  optional parameter to define a lower cutoff for the ALD length given in cM
# param: hval = optional parameter to define upper cutoff for the ALD length given in cM
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


input <- "Example_ALD.txt"

Output_prefix <- as.character("Example_ALD_output")  
lval <- as.numeric(0.05)	
hval <- as.numeric(1)		
constant <- T              
tm_lower <-  100
tm_upper <-  5000
k_lower <- 2
k_upper <- 1e8

Get_points <- function(input,lval,hval){
  # Read input file
  data <- read.table(input, header = F)
  
  
  # set dist and LD
  dist <- data[,1]
  LD <- data[,2]
  ndist <- length(dist)  ## number of rows in dataset
  lval=lval
  hval=hval
  
  # check x lower value and y lower value
  data.sub <- data
  if ((lval > dist[1]) || (hval < dist[ndist])) {
    data.sub <- subset(data, ((dist <= hval) & (dist >= lval)))
  }
  dist <- data.sub[,1]		# updated x values
  LD <- data.sub[,2]		# updated y values
  result_table <- data.frame(dist,LD)
  result_table <- result_table[complete.cases(result_table),]
  return(result_table)
}

Simple_Pulse_fn <- function(Data,constant,tm_lower,tm_upper){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$LD
  #Fitting an Exponential
  
  if(constant){
    fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])+x[3]
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,tm_lower,-1), upper=c(1, tm_upper,1), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "tm","c")
    
    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)+c), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE)) 
  } else {
    fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,tm_lower), upper=c(1, tm_upper), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "tm")
    
    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))
    
  }
  return(fit1_exp)
}

Get_RSS_fn <- function(SP_fit,EP_fit){
  SP_rss = sum(resid(SP_fit)^2)
  EP_rss = sum(resid(EP_fit)^2)
  return(data.frame(SP_rss=SP_rss,EP_rss=EP_rss))
  
}


Extended_Pulse_fn <- function(Data,constant,tm_lower,tm_upper,k_lower,k_upper){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$LD

  if(constant){
    fm1_lomax <- function(x) x[4] + x[3]* (1/(1 + (  (x[2] / x[1])*(dist/100) )))^(x[1])
    fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(k_lower,tm_lower,1e-6,-1), upper=c(k_upper,tm_upper,1,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],par1_lomax[2],par1_lomax[3],par1_lomax[4])
    names(par1_lomax) <- c("k","tm","A","c")
    fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((tm/ k) *(dist/100)) ))^(k), start=par1_lomax, algorithm="port",
                      lower=c(k_lower,tm_lower,1e-6,-1), upper=c(k_upper,tm_upper,1,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
    
  } else {
    fm1_lomax <- function(x)  x[3]* (1/(1 + (  (x[2] / x[1])*(dist/100) )))^(x[1])
    fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(k_lower,tm_lower,0), upper=c(k_upper,tm_upper,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],par1_lomax[2],par1_lomax[3])
    names(par1_lomax) <- c("k","tm","A")
    fit1_lomax <- nls(wcorr ~ A*(1/(1  + ((tm/ k) *(dist/100)) ))^(k), start=par1_lomax, algorithm="port",
                      lower=c(k_lower,tm_lower,0), upper=c(k_upper,tm_upper,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
    
    
    
  }
    

  
  return(fit1_lomax)
}


# Fitting the simple pulse

xdata <- Get_points(input,lval,hval)

SP_Fit <- Simple_Pulse_fn(xdata,constant,tm_lower,tm_upper)
A_SP_est <- as.numeric(coef(SP_Fit)[1])
tm_SP_est <- as.numeric(coef(SP_Fit)[2])  	# rate of decay of exponential
if(constant){
  C_SP_est <- as.numeric(coef(SP_Fit)[3])
  RSS <- sum(resid(SP_Fit)^2)
  
  Fit_simple_pulse_df <- data.frame(intercept = A_SP_est, tm = tm_SP_est,constant=C_SP_est, RSS = RSS)
  
}else{
  RSS <- sum(resid(SP_Fit)^2)
  
  Fit_simple_pulse_df <- data.frame(intercept = A_SP_est, tm = tm_SP_est, RSS = RSS)
}

write.csv(Fit_simple_pulse_df,file = paste(Output_prefix,"_simple_pulse.txt",sep = ""),quote = F,row.names = F)



# Fit the extended pulse

EP_Fit <- Extended_Pulse_fn(xdata,constant,tm_lower,tm_upper,k_lower,k_upper)
A_EP_est <- as.numeric(coef(EP_Fit)[3])
tm_EP_est <- as.numeric(coef(EP_Fit)[1])  	
k_EP_est <- as.numeric(coef(EP_Fit)[2])  
if(constant){
  C_EP_est <- as.numeric(coef(EP_Fit)[3])
  RSS <- sum(resid(EP_Fit)^2)
  
  Fit_extended_pulse_df <- data.frame(intercept = A_EP_est, tm = tm_EP_est,k = k_EP_est, constant= C_EP_est,RSS = RSS)
  
}else{
  RSS <- sum(resid(Fit_extended_pulse)^2)
  
  Fit_extended_pulse_df <- data.frame(intercept = A_EP_est, tm = tm_EP_est,k = k_EP_est,RSS = RSS)
  
}

write.csv(Fit_extended_pulse_df,file = paste(Output_prefix,"_extended_pulse.txt",sep = ""),quote = F,row.names = F)

