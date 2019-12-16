### Fit both Exponential and Lomax

suppressPackageStartupMessages({
  library(VGAM)
  library(tidyverse)
  library("DEoptim")
  library("MASS")
  library(bbmle)
})



input <- as.character(snakemake@input[[1]]) 		# output from rollof

Snake_output_One <- as.character(snakemake@output[[1]])  		# output of expfit log file
lval <- as.numeric(snakemake@wildcards[['min_dist_Fit']])		# lower value of dist to use
hval <- as.numeric(snakemake@params[['max_dist']])		# higher value of dist to use  

iterations=10
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

Fit_Exp_fn <- function(Data,affine){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$wcorr
  #Fitting an Exponential
  
  if(affine){
    fm1_exp <- function(x) x[1]*exp(-(dist/100)/x[2])+x[3]
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,0,-1), upper=c(1, 1,1), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "s","c")
    A_exp_est <- as.numeric((par1_exp[1]))
    Lambda_est <-as.numeric((par1_exp[2]))  	# rate of decay of exponential
    C_exp_est <- as.numeric((par1_exp[3]))
    
    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)/s)+c), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE)) 
  } else {
    fm1_exp <- function(x) x[1]*exp(-(dist/100)/x[2])
    fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
    fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,0), upper=c(1, 1), control=list(trace=FALSE))
    
    par1_exp <- fm3_exp$optim$bestmem
    # parameters for y ~ Ae-mt
    names(par1_exp) <- c("A", "s")
    A_exp_est <- as.numeric((par1_exp[1]))
    Lambda_est <-as.numeric((par1_exp[2]))  	# rate of decay of exponential
    
    fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)/s)), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE))
    
  }
  return(fit1_exp)
}



Fit_Lomax_fn <- function(Data,affine,Expo_s){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$wcorr
  if(affine){
    fm1_lomax <- function(x) x[4] + x[3]* (1/(1 + (x[1]*(dist/100) /  x[2])))^(1/x[1])
    fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1e-6,Expo_s,0,0), upper=c(100,Expo_s,1,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],Expo_s,par1_lomax[3],par1_lomax[4])
    names(par1_lomax) <- c("w","s","A","c")
    fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((w*(dist/100)) / s)))^(1/w), start=par1_lomax, control=list(maxiter=10000, warnOnly=TRUE))
    
  } else {
    fm1_lomax <- function(x)  x[3]* (1/(1 + (x[1]*(dist/100) /  x[2])))^(1/x[1])
    fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
    fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(1e-6,Expo_s/2,0), upper=c(100,Expo_s*2,1), control=list(trace=FALSE))
    
    par1_lomax <- fm3_DEoptim$optim$bestmem
    par1_lomax <- c(par1_lomax[1],Expo_s,par1_lomax[3])
    names(par1_lomax) <- c("w","s","A")
    fit1_lomax <- nls(wcorr ~ A*(1/(1  + ((w*(dist/100)) / s)))^(1/w), start=par1_lomax, control=list(maxiter=10000, warnOnly=TRUE))
    
  }
  return(fit1_lomax)
}

Fit_Lomax_multiple_times <- function(iterations,xdata,affine=T,Expo_s){
  list_Fits <- list()
  Fit_name <- c()
  for(i in 1:iterations){
    alternativeFunction <- function(xx){return(list(NA))}
    options(warn = 1)
    End_Fit<-  try(Fit_Lomax_fn(Data = xdata,affine = affine,Expo_s=Expo_s),silent = FALSE)
    if (inherits(End_Fit, "try-error")) End_Fit=alternativeFunction(xdata)
    list_Fits[[i]] <- End_Fit
    print(paste('Fit iteration :',i,sep = " "))
    Fit_name_x <- c(paste("Fit_",i,sep=""))
    Fit_name <- c(Fit_name,Fit_name_x)
  }
  
  # set names
  names(list_Fits) = Fit_name
  
  
  # Calculat RSS for every fit
  RSS_comparison <- c()
  for(i in 1:iterations){
    if(is.na(list_Fits[[i]])==T){
      RSS <- 10
    } 
    else{
      RSS <- sum(resid(list_Fits[[i]])^2)
    }
    
    xx <- cbind(paste("Fit_",i,sep=""),as.numeric(as.character(RSS)))
    RSS_comparison <- rbind(RSS_comparison,xx)
  }
  
  RSS_comparison <- as.data.frame(RSS_comparison)
  RSS_comparison$V2 <- as.numeric(as.character(RSS_comparison$V2))
  RSS_comparison$V1 <- as.character(RSS_comparison$V1)
  #Choose the fit with the smallest RSS
  xx=RSS_comparison$V1[RSS_comparison$V2==min(RSS_comparison$V2)]
  Best_fitting_Lomax_model <- list_Fits[[xx[1]]]
  return(Best_fitting_Lomax_model)
}

xdata <- Get_points(input,lval,hval,log = F)

Fit_Exp <- Fit_Exp_fn(xdata,affine = T)

s_exp_est <- as.numeric(coef(Fit_Exp)[2]) 




options(warn = 1)
Best_fitting_Lomax_model <-  try(Fit_Lomax_multiple_times(iterations,xdata,affine=T,Expo_s=s_exp_est),silent = FALSE)
if (inherits(Best_fitting_Lomax_model, "try-error")){
  if(hval == 1) {
    xdata <- Get_points(input,lval,hval+1,log = F)
    
    Fit_Exp <- Fit_Exp_fn(xdata,affine = T)
    
    s_exp_est <- as.numeric(coef(Fit_Exp)[2]) 
    
    options(warn = 1)
    Best_fitting_Lomax_model_varying_hval <-  try(Fit_Lomax_multiple_times(iterations,xdata,affine=T,Expo_s=s_exp_est),silent = FALSE)
    if (inherits(Best_fitting_Lomax_model_varying_hval, "try-error")){
      A_exp_est <- as.numeric(coef(Fit_Exp)[1])
      s_exp_est <- as.numeric(coef(Fit_Exp)[2])  	# rate of decay of exponential
      C_exp_est <- as.numeric(coef(Fit_Exp)[3])
      Test_Sig=NA
      
      #print output
      outlog <- paste(Snake_output_One)
      cat("Summary of fit:\n", file=outlog, append = FALSE)
      
      # Convert distance in Morgans
      capture.output(Test_Sig, file = outlog, append = TRUE)
      cat(paste("A, s, c, RSS_Expo, AIC_Expo, A, s, w,c, RSS_Lomax, AIC_Lomax, F_Test:","\t",A_exp_est,"\t",s_exp_est,"\t",C_exp_est,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\n",sep = ""), file=outlog, append=TRUE)
      
    } else{
      
      Test_Sig=anova(Fit_Exp,Best_fitting_Lomax_model_varying_hval,test="F")
      AIC_Test <- AICtab(Fit_Exp,Best_fitting_Lomax_model_varying_hval)
      select.exp=attr(AIC_Test,'row.names')=="Fit_Exp"
      if (select.exp[1]==T){
        Expo.AIC=AIC_Test$dAIC[1]
        Lomax.AIC=AIC_Test$dAIC[2]
      }else {
        Expo.AIC=AIC_Test$dAIC[2]
        Lomax.AIC=AIC_Test$dAIC[1]
      }
      
      
      A_exp_est <- as.numeric(coef(Fit_Exp)[1])
      s_exp_est <- as.numeric(coef(Fit_Exp)[2])  	# rate of decay of exponential
      C_exp_est <- as.numeric(coef(Fit_Exp)[3])
      
      A_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[3])
      s_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[2])  	# rate of decay of exponential
      w_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[1])
      C_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[4])
      
      #print output
      outlog <- paste(Snake_output_One)
      cat("Summary of fit:\n", file=outlog, append = FALSE)
      
      # Convert distance in Morgans
      capture.output(Test_Sig, file = outlog, append = TRUE)
      cat(paste("A, s, c, RSS_Expo, AIC_Expo, A, s, w,c, RSS_Lomax, AIC_Lomax, F_Test:","\t",A_exp_est,"\t",s_exp_est,"\t",C_exp_est,"\t",Test_Sig$`Res.Sum Sq`[1],"\t",Expo.AIC,"\t",A_lomax_est,"\t",s_lomax_est,"\t",w_lomax_est,"\t",C_lomax_est,"\t",Test_Sig$`Res.Sum Sq`[2],"\t",Lomax.AIC,"\t",Test_Sig$`Pr(>F)`[2],"\n",sep = ""), file=outlog, append=TRUE)
      
    }
  }
  if(hval > 1) {
    
    Fit_all_using_diff_hval <- function(iterations,input,lval,hval,log,affine=T){
      
      xdata <- Get_points(input,lval,hval+1,log = F)
      Fit_Exp <- Fit_Exp_fn(xdata,affine = T)
      s_exp_est <- as.numeric(coef(Fit_Exp)[2])
      print(s_exp_est)
      list_Fits <- list()
      Fit_name <- c()
      for(i in 1:(iterations/2)){
        alternativeFunction <- function(xx){return(list(NA))}
        options(warn = 1)
        End_Fit<-  try(Fit_Lomax_fn(Data = xdata,affine = affine,Expo_s=s_exp_est),silent = FALSE)
        if (inherits(End_Fit, "try-error")) End_Fit=alternativeFunction(xdata)
        list_Fits[[i]] <- End_Fit
        print(paste('Fit iteration :',i,sep = " "))
        Fit_name_x <- c(paste("Fit_",i,sep=""))
        Fit_name <- c(Fit_name,Fit_name_x)
      }
      
      xdata <- Get_points(input,lval,hval-1,log = F)
      Fit_Exp <- Fit_Exp_fn(xdata,affine = T)
      s_exp_est <- as.numeric(coef(Fit_Exp)[2])
      
      for(i in 1:(iterations/2)){
        alternativeFunction <- function(xx){return(list(NA))}
        options(warn = 1)
        End_Fit<-  try(Fit_Lomax_fn(Data = xdata,affine = affine,Expo_s=s_exp_est),silent = FALSE)
        if (inherits(End_Fit, "try-error")) End_Fit=alternativeFunction(xdata)
        list_Fits[[(i+(iterations/2))]] <- End_Fit
        print(paste('Fit iteration :',i,sep = " "))
        Fit_name_x <- c(paste("Fit_",i+(iterations/2),sep=""))
        Fit_name <- c(Fit_name,Fit_name_x)
      }
      
      # set names
      names(list_Fits) = Fit_name
      
      
      # Calculat RSS for every fit
      RSS_comparison <- c()
      for(i in 1:iterations){
        if(is.na(list_Fits[[i]])==T){
          RSS <- 10
        } 
        else{
          RSS <- sum(resid(list_Fits[[i]])^2)
        }
        
        xx <- cbind(paste("Fit_",i,sep=""),as.numeric(as.character(RSS)))
        RSS_comparison <- rbind(RSS_comparison,xx)
      }
      
      RSS_comparison <- as.data.frame(RSS_comparison)
      RSS_comparison$V2 <- as.numeric(as.character(RSS_comparison$V2))
      RSS_comparison$V1 <- as.character(RSS_comparison$V1)
      #Choose the fit with the smallest RSS
      xx=RSS_comparison$V1[RSS_comparison$V2==min(RSS_comparison$V2)]
      Best_fitting_Lomax_model <- list_Fits[[xx[1]]]
      print(RSS_comparison)
      return(Best_fitting_Lomax_model)
    }
    options(warn = 1)
    Best_fitting_Lomax_model_varying_hval <-  try(Fit_all_using_diff_hval(iterations,input,lval,hval,log=F,affine=T),silent = FALSE)
    if (inherits(Best_fitting_Lomax_model_varying_hval, "try-error")){
      A_exp_est <- as.numeric(coef(Fit_Exp)[1])
      s_exp_est <- as.numeric(coef(Fit_Exp)[2])  	# rate of decay of exponential
      C_exp_est <- as.numeric(coef(Fit_Exp)[3])
      Test_Sig=NA
      
      #print output
      outlog <- paste(Snake_output_One)
      cat("Summary of fit:\n", file=outlog, append = FALSE)
      
      # Convert distance in Morgans
      capture.output(Test_Sig, file = outlog, append = TRUE)
      cat(paste("A, s, c, RSS_Expo, AIC_Expo, A, s, w,c, RSS_Lomax, AIC_Lomax, F_Test:","\t",A_exp_est,"\t",s_exp_est,"\t",C_exp_est,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\t",NA,"\n",sep = ""), file=outlog, append=TRUE)
      
    } else{
      
      Test_Sig=anova(Fit_Exp,Best_fitting_Lomax_model_varying_hval,test="F")
      AIC_Test <- AICtab(Fit_Exp,Best_fitting_Lomax_model_varying_hval)
      select.exp=attr(AIC_Test,'row.names')=="Fit_Exp"
      if (select.exp[1]==T){
        Expo.AIC=AIC_Test$dAIC[1]
        Lomax.AIC=AIC_Test$dAIC[2]
      }else {
        Expo.AIC=AIC_Test$dAIC[2]
        Lomax.AIC=AIC_Test$dAIC[1]
      }
      
      
      A_exp_est <- as.numeric(coef(Fit_Exp)[1])
      s_exp_est <- as.numeric(coef(Fit_Exp)[2])  	# rate of decay of exponential
      C_exp_est <- as.numeric(coef(Fit_Exp)[3])
      
      A_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[3])
      s_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[2])  	# rate of decay of exponential
      w_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[1])
      C_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model_varying_hval)[4])
      
      #print output
      outlog <- paste(Snake_output_One)
      cat("Summary of fit:\n", file=outlog, append = FALSE)
      
      # Convert distance in Morgans
      capture.output(Test_Sig, file = outlog, append = TRUE)
      cat(paste("A, s, c, RSS_Expo, AIC_Expo, A, s, w,c, RSS_Lomax, AIC_Lomax, F_Test:","\t",A_exp_est,"\t",s_exp_est,"\t",C_exp_est,"\t",Test_Sig$`Res.Sum Sq`[1],"\t",Expo.AIC,"\t",A_lomax_est,"\t",s_lomax_est,"\t",w_lomax_est,"\t",C_lomax_est,"\t",Test_Sig$`Res.Sum Sq`[2],"\t",Lomax.AIC,"\t",Test_Sig$`Pr(>F)`[2],"\n",sep = ""), file=outlog, append=TRUE)
    }
    
  }
    
  } else{
    
    Test_Sig=anova(Fit_Exp,Best_fitting_Lomax_model,test="F")
    AIC_Test <- AICtab(Fit_Exp,Best_fitting_Lomax_model)
    select.exp=attr(AIC_Test,'row.names')=="Fit_Exp"
    if (select.exp[1]==T){
      Expo.AIC=AIC_Test$dAIC[1]
      Lomax.AIC=AIC_Test$dAIC[2]
    }else {
      Expo.AIC=AIC_Test$dAIC[2]
      Lomax.AIC=AIC_Test$dAIC[1]
    }
    
    
    A_exp_est <- as.numeric(coef(Fit_Exp)[1])
    s_exp_est <- as.numeric(coef(Fit_Exp)[2])  	# rate of decay of exponential
    C_exp_est <- as.numeric(coef(Fit_Exp)[3])
    
    A_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model)[3])
    s_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model)[2])  	# rate of decay of exponential
    w_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model)[1])
    C_lomax_est <- as.numeric(coef(Best_fitting_Lomax_model)[4])
    
    #print output
    outlog <- paste(Snake_output_One)
    cat("Summary of fit:\n", file=outlog, append = FALSE)
    
    # Convert distance in Morgans
    capture.output(Test_Sig, file = outlog, append = TRUE)
    cat(paste("A, s, c, RSS_Expo, AIC_Expo, A, s, w,c, RSS_Lomax, AIC_Lomax, F_Test:","\t",A_exp_est,"\t",s_exp_est,"\t",C_exp_est,"\t",Test_Sig$`Res.Sum Sq`[1],"\t",Expo.AIC,"\t",A_lomax_est,"\t",s_lomax_est,"\t",w_lomax_est,"\t",C_lomax_est,"\t",Test_Sig$`Res.Sum Sq`[2],"\t",Lomax.AIC,"\t",Test_Sig$`Pr(>F)`[2],"\n",sep = ""), file=outlog, append=TRUE)
    
  }
  