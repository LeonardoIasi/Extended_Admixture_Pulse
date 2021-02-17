### Fit the simple and extended pulse model to ALD data

## The ALD data should be in the raw ALDER output format as follows:
# col 1 = distance between bins of SNPs the LD is computed for in centiMorgan
# col 2 = the LD per bin
# col 3 (optional) = bin count

suppressPackageStartupMessages({
  library(VGAM)
  library("DEoptim")
  library("MASS")

})



input <- as.character(snakemake@input[[1]]) 		# output from rollof/ALDER 

Output <- as.character("test_output.txt")  		# output file name
lval <- as.numeric(0.05)		# lower value of distance to use
hval <- as.numeric(10)		# higher value of dist to use
affine <- T               # sets c as parameter

Get_points <- function(input,lval,hval,log){
  # Read input file
  data <- read.table(input, header = F)
  
  
  # set dist and LD
  col=2
  dist <- data[,1]
  LD <- data[,col]
  ndist <- length(dist)  ## number of rows in dataset
  lval=lval
  hval=hval
  
  # check x lower value and y lower value
  data.sub <- data
  if ((lval > dist[1]) || (hval < dist[ndist])) {
    data.sub <- subset(data, ((dist <= hval) & (dist >= lval)))
  }
  dist <- data.sub[,1]		# updated x values
  LD <- data.sub[,col]		# updated y values
  if(log==T){
    xx <- cbind(dist,log(LD))
    xx <- xx[complete.cases(xx),]
    LD <- xx[,2]
    LD <-c(LD,rep(NA,length(data.sub[,1])-length(LD)))
    dist <- xx[,1]
    dist <-c(dist,rep(NA,length(data.sub[,1])-length(dist)))
  }
  result_table <- data.frame(dist,LD)
  result_table <- result_table[complete.cases(result_table),]
  return(result_table)
}

Fit_simple_p_fn <- function(Data,affine){
  xx=Data
  
  dist=xx$dist
  LD=xx$LD
  #Fitting an Exponential
  
  if(affine){
    fm1_exp <- function(x) x[1]*exp(-(dist/100)/x[2])+x[3]
    fm2_exp <- function(x) sum((LD-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,0,-1), upper=c(1, 1,1), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "lambda","c")
    A_exp_est <- as.numeric((par1_exp[1]))
    Lambda_est <-as.numeric((par1_exp[2]))  	# rate of decay of exponential
    C_exp_est <- as.numeric((par1_exp[3]))
    
    fit1_exp <- nls(LD ~ (A*exp(-(dist/100)/lambda)+c), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))
  } else {
    fm1_exp <- function(x) x[1]*exp(-(dist/100)/x[2])
    fm2_exp <- function(x) sum((LD-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,0), upper=c(1, 1), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "lambda")
    A_exp_est <- as.numeric((par1_exp[1]))
    Lambda_est <-as.numeric((par1_exp[2]))  	# rate of decay of exponential
    
    fit1_exp <- nls(LD ~ (A*exp(-(dist/100)/lambda)), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))
    
  }
  return(fit1_exp)
}

### uses lomax with shape = alpha (k) and scale = lambda (tm/k)

Fit_extended_pulse_fn <- function(Data,affine){
  xx=Data
  
  dist=xx$dist
  LD=xx$LD
  # here lambda = (tm/k) and alpha = k
  if(affine){
    fm1_lomax <- function(x) x[4] + x[3]* (1 + (dist/100) / x[1])^(-x[2])
    fm2_lomax_k <- function(x) sum((LD-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1e-6,1,0,0), upper=c(10,100,1,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],par1_lomax[2],par1_lomax[3],par1_lomax[4])
    names(par1_lomax) <- c("lambda","alpha","A","c")
    fit1_lomax <- nls(LD ~ c+A*(1 + (dist/100)/ lambda )^(-alpha), start=par1_lomax, control=list(maxiter=10000, warnOnly=TRUE,minFactor=0.0004))
    
  } else {
    fm1_lomax <- function(x) x[3]* (1 + (dist/100) / x[1])^(-x[2])
    fm2_lomax_k <- function(x) sum((LD-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1e-6,Expo_s/2,0), upper=c(100,Expo_s*2,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],par1_lomax[2],par1_lomax[3])
    names(par1_lomax) <- c("lambda","alpha","A")
    fit1_lomax <- nls(LD ~  A*(1 + (dist/100)/ lambda )^(-alpha), start=par1_lomax, control=list(maxiter=10000, warnOnly=TRUE,minFactor=0.0004))
    
  }
  
  
  
  return(fit1_lomax)
}


# Fitting the simple pulse

xdata <- Get_points(input,lval,hval,log = F)

Fit_simple_pulse <- Fit_simple_p_fn(xdata,affine = T)
A_SP_est <- as.numeric(coef(Fit_simple_pulse)[1])
Lambda_SP_est <- as.numeric(coef(Fit_simple_pulse)[2])  	# rate of decay of exponential
if(affline){
  C_SP_est <- as.numeric(coef(Fit_simple_pulse)[3])
  RSS <- sum(resid(Fit_simple_pulse)^2)
  
  Fit_simple_pulse_df <- data.frame(intercept = A_SP_est, lambda = Lambda_SP_est,affline=C_SP_est, RSS = RSS)
  
}else{
  RSS <- sum(resid(Fit_simple_pulse)^2)
  
  Fit_simple_pulse_df <- data.frame(intercept = A_SP_est, lambda = Lambda_SP_est, RSS = RSS)
  
}

write.csv(Fit_simple_pulse_df)



# Fit the extended pulse

Fit_extended_pulse <- Fit_extended_pulse_fn(xdata,affine = T)
A_EP_est <- as.numeric(coef(Fit_extended_pulse)[3])
Lambda_EP_est <- as.numeric(coef(Fit_extended_pulse)[1])  	
alpha_EP_est <- as.numeric(coef(Fit_extended_pulse)[2])  
if(affline){
  C_EP_est <- as.numeric(coef(Fit_extended_pulse)[3])
  RSS <- sum(resid(Fit_Exp)^2)
  
  Fit_extended_pulse_df <- data.frame(intercept = A_EP_est, lambda = Lambda_EP_est,alpha = alpha_EP_est, affline= C_EP_est,RSS = RSS)
  
}else{
  RSS <- sum(resid(Fit_extended_pulse)^2)
  
  Fit_extended_pulse_df <- data.frame(intercept = A_EP_est, lambda = Lambda_EP_est,alpha = alpha_EP_est, RSS = RSS)
  
}

write.csv(Fit_extended_pulse)

