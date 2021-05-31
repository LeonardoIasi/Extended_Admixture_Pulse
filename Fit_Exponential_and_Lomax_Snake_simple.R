### Fit both Exponential and Lomax

suppressPackageStartupMessages({
  library(VGAM)
  library("DEoptim")
  library("MASS")
  library(bbmle)
})


input <- as.character(snakemake@input[[1]]) 		# output from rollof

Snake_output_One <- as.character(snakemake@output[[1]])  		# output of expfit log file
lval <- as.numeric(snakemake@wildcards[['min_dist_Fit']])		# lower value of dist to use
hval <- as.numeric(snakemake@params[['max_dist']])		# higher value of dist to use
Fix_tm <- as.logical(snakemake@params[['Fix_lambda']])		# higher value of dist to use
Only_Simple_Pulse_fit <- as.logical(snakemake@params[['Only_Simple_Pulse_fit']]) # if only a simple pulse should be fitted

Get_points <- function(input,lval,hval,log){
  # Read input file
  data <- read.table(input, header = F)


  # set dist and wcorr
  col=2
  dist <- data[,1]
  wcorr <- data[,col]
  ndist <- length(dist)  ## number of rows in dataset
  lval=lval
  hval=hval

  # check x lower value and y lower value
  data.sub <- data
  if ((lval > dist[1]) || (hval < dist[ndist])) {
    data.sub <- subset(data, ((dist <= hval) & (dist >= lval)))
  }
  dist <- data.sub[,1]		# updated x values
  wcorr <- data.sub[,col]		# updated y values
  if(log==T){
    xx <- cbind(dist,log(wcorr))
    xx <- xx[complete.cases(xx),]
    wcorr <- xx[,2]
    wcorr <-c(wcorr,rep(NA,length(data.sub[,1])-length(wcorr)))
    dist <- xx[,1]
    dist <-c(dist,rep(NA,length(data.sub[,1])-length(dist)))
  }
  result_table <- data.frame(dist,wcorr)
  result_table <- result_table[complete.cases(result_table),]
  return(result_table)
}

Simple_Pulse_fn <- function(Data,affine){
  xx=Data

  dist=xx$dist
  wcorr=xx$wcorr
  #Fitting an Exponential

  if(affine){
    fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])+x[3]
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,1,-1), upper=c(1, 5000,1), control=list(trace=FALSE))

    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "tm","c")
    A_exp_est <- as.numeric((par1_exp[1]))
    Lambda_est <-as.numeric((par1_exp[2]))  	# rate of decay of exponential
    C_exp_est <- as.numeric((par1_exp[3]))

    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)+c), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))
  } else {
    fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,1), upper=c(1, 5000), control=list(trace=FALSE))

    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "tm")
    A_exp_est <- as.numeric((par1_exp[1]))
    Lambda_est <-as.numeric((par1_exp[2]))  	# rate of decay of exponential

    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))

  }
  return(fit1_exp)
}

Get_RSS_fn <- function(SP_fit,EP_fit){
  SP_rss = sum(resid(SP_fit)^2)
  EP_rss = sum(resid(EP_fit)^2)
  return(data.frame(SP_rss=SP_rss,EP_rss=EP_rss))

}



AIC_Test_fn <- function(Fit_SP,Fit_EP){
  AIC_Test = bbmle::AICtab(Fit_SP,Fit_EP)
  select.exp=attr(AIC_Test,'row.names')=="Fit_SP"
  if (select.exp[1]==T){
    SP.AIC=AIC_Test$dAIC[1]
    EP.AIC=AIC_Test$dAIC[2]
  }else {
    SP.AIC=AIC_Test$dAIC[2]
    EP.AIC=AIC_Test$dAIC[1]
  }
  return(data.frame(SP.AIC=SP.AIC,EP.AIC=EP.AIC))
}

Extended_Pulse_fn <- function(Data,affine,Expo_tm,Fix_tm){
  xx=Data

  dist=xx$dist
  wcorr=xx$wcorr
  if(Fix_tm){
    print("Using fixed mean time")
    if(affine){
      fm1_lomax <- function(x) x[4] + x[3]* (1/(1 + (  (x[2] / x[1])*(dist/100) )))^(x[1])
      fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
      fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1,Expo_tm,0,0), upper=c(1e8,Expo_tm,1,1), control=list(trace=FALSE))

      par1_lomax <- fm3_DEoptim$optim$bestmem
      par1_lomax <- c(par1_lomax[1],par1_lomax[3],par1_lomax[4])
      names(par1_lomax) <- c("k","A","c")
      fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((Expo_tm/ k) *(dist/100)) ))^(k), start=par1_lomax, algorithm="port",
                        lower=c(1,0,-1), upper=c(1e8,1,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))

    } else {
      fm1_lomax <- function(x) x[3]* (1/(1 + (  (x[2] / x[1])*(dist/100) )))^(x[1])
      fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
      fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1,Expo_tm,0), upper=c(10000,Expo_tm,1), control=list(trace=FALSE))

      par1_lomax <- fm3_DEoptim$optim$bestmem
      par1_lomax <- c(par1_lomax[1],par1_lomax[3])
      names(par1_lomax) <- c("k","A")
      fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((Expo_tm/ k) *(dist/100)) ))^(k), start=par1_lomax, algorithm="port",
                        lower=c(1,0), upper=c(1e8,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
      #fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((Expo_tm/ k) *(dist/100)) ))^(k), start=par1_lomax, control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))

    }
  } else{
    print("Estimating tm")
    if(affine){
      fm1_lomax <- function(x) x[4] + x[3]* (1/(1 + (  (x[2] / x[1])*(dist/100) )))^(x[1])
      fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
      fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(2,100,1e-6,-1), upper=c(1e8,5000,1,1), control=list(trace=FALSE))

      par1_lomax <- fm3_DEoptim$optim$bestmem
      par1_lomax <- c(par1_lomax[1],Expo_tm,par1_lomax[3],par1_lomax[4])
      names(par1_lomax) <- c("k","tm","A","c")
      fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((tm/ k) *(dist/100)) ))^(k), start=par1_lomax, algorithm="port",
                        lower=c(2,100,0,-1), upper=c(1e8,5000,1,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
      # use Gauss-Newton algorithm
      #fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((tm/ k) *(dist/100)) ))^(k), start=par1_lomax, control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))

    } else {
            fm1_lomax <- function(x)  x[3]* (1/(1 + (  (x[2] / x[1])*(dist/100) )))^(x[1])
      fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
      fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(2,100,0), upper=c(1e8,5000,1), control=list(trace=FALSE))

      par1_lomax <- fm3_DEoptim$optim$bestmem
      par1_lomax <- c(par1_lomax[1],Expo_tm,par1_lomax[3])
      names(par1_lomax) <- c("k","tm","A")
      fit1_lomax <- nls(wcorr ~ A*(1/(1  + ((tm/ k) *(dist/100)) ))^(k), start=par1_lomax, algorithm="port",
                        lower=c(2,100,0), upper=c(1e8,5000,1),control=list(maxiter=100000, warnOnly=TRUE,minFactor=0.0004))
      # use Gauss-Newton algorithm
      fit1_lomax <- nls(wcorr ~ A*(1/(1  + ((tm/ k) *(dist/100)) ))^(k), start=par1_lomax, control=list(maxiter=10000, warnOnly=TRUE))
      


    }

  }

  return(fit1_lomax)
}

xdata <- Get_points(input,lval,hval,log = F)

options(warn = 1)
Fit_SP <-  try(Simple_Pulse_fn(xdata,affine = T),silent = FALSE)
if (inherits(Fit_SP, "try-error")){
  print("Fitting failed")
  
  #print output
  outlog <- paste(Snake_output_One)
  cat("Summary of fit:\n", file=outlog, append = FALSE)
  F_Test=NA
  # Convert distance in Morgans
  capture.output(F_Test, file = outlog, append = TRUE)
  cat(paste("A, tm, c, RSS_SP, AIC_SP, A, tm, k,c, RSS_EP, AIC_EP, F_Test:","\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\n",sep = ""), file=outlog, append=TRUE)
}  else{
  Expo_tm <- as.numeric(coef(Fit_SP)[2])
  if(Only_Simple_Pulse_fit){
    A_sp_est <- as.numeric(coef(Fit_SP)[1])
    tm_sp_est <- as.numeric(coef(Fit_SP)[2])  	# rate of decay of exponential
    C_sp_est <- as.numeric(coef(Fit_SP)[3])
    F_Test=NA
    
    #print output
    outlog <- paste(Snake_output_One)
    cat("Summary of fit:\n", file=outlog, append = FALSE)
    
    # Convert distance in Morgans
    capture.output(F_Test, file = outlog, append = TRUE)
    cat(paste("A, tm, c, RSS_SP, AIC_SP, A, tm, k,c, RSS_EP, AIC_EP, F_Test:","\t",A_sp_est,"\t",tm_sp_est,"\t",C_sp_est,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\n",sep = ""), file=outlog, append=TRUE)
    
  }else{
    options(warn = 1)
    Fit_EP <-  try(Extended_Pulse_fn(xdata,affine=T,Expo_tm=Expo_tm,Fix_tm),silent = FALSE)
    if (inherits(Fit_EP, "try-error")){
      
      A_sp_est <- as.numeric(coef(Fit_SP)[1])
      tm_sp_est <- as.numeric(coef(Fit_SP)[2])  	# rate of decay of exponential
      C_sp_est <- as.numeric(coef(Fit_SP)[3])
      F_Test=NA
      
      #print output
      outlog <- paste(Snake_output_One)
      cat("Summary of fit:\n", file=outlog, append = FALSE)
      
      # Convert distance in Morgans
      capture.output(F_Test, file = outlog, append = TRUE)
      cat(paste("A, tm, c, RSS_SP, AIC_SP, A, tm, k,c, RSS_EP, AIC_EP, F_Test:","\t",A_sp_est,"\t",tm_sp_est,"\t",C_sp_est,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\n",sep = ""), file=outlog, append=TRUE)
      
    } else{
      
      A_sp_est <- as.numeric(coef(Fit_SP)[1])
      tm_sp_est <- as.numeric(coef(Fit_SP)[2])  	# rate of decay of exponential
      C_sp_est <- as.numeric(coef(Fit_SP)[3])
      
      A_ep_est <- as.numeric(coef(Fit_EP)[3])
      if(Fix_tm){tm_ep_est <- as.numeric(coef(Fit_SP)[2])}
      else{tm_ep_est <- as.numeric(coef(Fit_EP)[2])}
      k_ep_est <- as.numeric(coef(Fit_EP)[1])
      C_ep_est <- as.numeric(coef(Fit_EP)[4])
      
      Compare_RSS = Get_RSS_fn(SP_fit = Fit_SP,EP_fit = Fit_EP)
      F_Test = anova(Fit_SP,Fit_EP,test="F")
      Compare_AIC = AIC_Test_fn(Fit_SP = Fit_SP,Fit_EP = Fit_EP)
      
      #print output
      outlog <- paste(Snake_output_One)
      cat("Summary of fit:\n", file=outlog, append = FALSE)
      
      # Convert distance in Morgans
      capture.output(F_Test, file = outlog, append = TRUE)
      cat(paste("A, tm, c, RSS_SP, AIC_SP, A, tm, k,c, RSS_EP, AIC_EP, F_Test:","\t",A_sp_est,"\t",tm_sp_est,"\t",C_sp_est,"\t",Compare_RSS$SP_rss,"\t",Compare_AIC$SP.AIC,"\t",A_ep_est,"\t",tm_ep_est,"\t",k_ep_est,"\t",C_ep_est,"\t",Compare_RSS$EP_rss,"\t",Compare_AIC$EP.AIC,"\t",F_Test$`Pr(>F)`[2],"\n",sep = ""), file=outlog, append=TRUE)
      
      
    }
  }
}




