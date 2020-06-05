# 1kG YRI CUE and the 3 hc Neandertals ALDER computed ALD Analysis

suppressPackageStartupMessages({
  library(VGAM)
  library(tidyverse)
  library("DEoptim")
  library("MASS")
  library(bbmle)
  library(viridis)
})

setwd("~/Desktop/ATE_Paper/Admixture_Time_Inference_Paper/Real_Data_Analysis/")
input <- "Raw_ALDER_output_ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_rr_CEU_specific_map.txt"
#input <- "Raw_ALDER_output-ALL_chr_hcNea_YRI_CEU_1k_G_LES_ascertained_corrected_constant_rr.txt"
lval <- 0.05	# lower value of dist to use
hval <- 3		# higher value of dist to use 
log <- F
affine <- T

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

Fit_Lomax_fn <- function(Data,fixed_w){
  xx=Data
  
  dist=xx$dist
  wcorr=xx$wcorr

  fm1_lomax <- function(x) x[4] + x[3]* (1/(1 + (x[1]*(dist/100) /  x[2])))^(1/x[1])
  fm2_lomax_k <- function(x) sum((wcorr-fm1_lomax(x))^2)
  fm3_DEoptim <- DEoptim(fm2_lomax_k, lower=c(fixed_w,0,0,0), upper=c(fixed_w,1,1,1), control=list(trace=FALSE))
  
  par1_lomax <- fm3_DEoptim$optim$bestmem
  par1_lomax <- c(par1_lomax[2],par1_lomax[3],par1_lomax[4])
  names(par1_lomax) <- c("s","A","c")
  fit1_lomax <- nls(wcorr ~ c+A*(1/(1  + ((fixed_w*(dist/100)) / s)))^(1/fixed_w), start=par1_lomax, control=list(maxiter=10000, warnOnly=TRUE,minFactor=0.0004))

  return(fit1_lomax)
}

xdata=Get_points(input,lval,hval,log)
Expo_fit_result <- Fit_Exp_fn(Data = xdata,affine = affine)
t_expo_est=1/coef(Expo_fit_result)[[2]]


t_d=c(1,100,200,400,800,1000,1500,2000,2500)

all_t_lomax_est <- list()
for(t in 1:length(t_d)){
  fixed_w= 1/(t_expo_est*(t_expo_est/((t_d[t]/4)^2)))
  
  Lomax_fit_result <- Fit_Lomax_fn(Data = xdata,fixed_w = fixed_w)
  coef(Lomax_fit_result)
  t_lomax_est=1/coef(Lomax_fit_result)[[1]]
  all_t_lomax_est[[t]] <- Lomax_fit_result
  
}

Mean_time_est <- data.frame(Model=c("Exp","Lomax(td=1)","Lomax(td=100)","Lomax(td=200)","Lomax(td=400)"
                                    ,"Lomax(td=800)","Lomax(td=1000)","Lomax(td=1500)","Lomax(td=2000)","Lomax(td=2500)"),
                            Mean_time=c(1/coef(Expo_fit_result)[[2]],1/coef(all_t_lomax_est[[1]])[[1]],1/coef(all_t_lomax_est[[2]])[[1]]
                                        ,1/coef(all_t_lomax_est[[3]])[[1]],1/coef(all_t_lomax_est[[4]])[[1]],1/coef(all_t_lomax_est[[5]])[[1]]
                                        ,1/coef(all_t_lomax_est[[6]])[[1]],1/coef(all_t_lomax_est[[7]])[[1]],1/coef(all_t_lomax_est[[8]])[[1]]
                                        ,1/coef(all_t_lomax_est[[9]])[[1]]))


EXP_Log_fn <- function(dist,A,s) log(A)  -((dist/100)/s)
LOMAX_Log_fn <- function(dist,A,s,w) log(A) + (1/w) * -log(1+ ((w*(dist/100))/s))


AIC_result <- AIC(Expo_fit_result,all_t_lomax_est[[1]],all_t_lomax_est[[2]],all_t_lomax_est[[3]],all_t_lomax_est[[4]],
    all_t_lomax_est[[5]],all_t_lomax_est[[6]],all_t_lomax_est[[7]],all_t_lomax_est[[8]],all_t_lomax_est[[9]])


Model_RSS <- c()
for(i in list(Expo_fit_result,all_t_lomax_est[[1]],all_t_lomax_est[[2]],all_t_lomax_est[[3]],all_t_lomax_est[[4]],
           all_t_lomax_est[[5]],all_t_lomax_est[[6]],all_t_lomax_est[[7]],all_t_lomax_est[[8]],all_t_lomax_est[[9]])){
  RSS <- sum(residuals(i)^2) 
  Model_RSS <- c(Model_RSS,RSS)
}

pdf("Progress_figures.pdf")
plot(xdata$dist,log(xdata$wcorr),pch=19,xlim=c(0,0.5))
lines(xdata$dist,EXP_Log_fn(xdata$dist,coef(Expo_fit_result)[[1]],coef(Expo_fit_result)[[2]]),col="red",lwd=2)
for(t in 1:length(t_d)){
  lines(xdata$dist,LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[t]])[[2]],coef(all_t_lomax_est[[t]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[t]/4)^2)))),col="green", lwd=2)
}
barplot(Mean_time_est$Mean_time,names.arg = Mean_time_est$Model,las=2,main="mean time est")
barplot(Model_RSS,names.arg = Mean_time_est$Model,las=2,main="RSS")
barplot(AIC_result$AIC+44000,names.arg = Mean_time_est$Model,las=2,main="AIC")
dev.off()

log_xdata <- data.frame(dist=xdata$dist,log_wcorr=log(xdata$wcorr))
log_xdata <- log_xdata[complete.cases(log_xdata$log_wcorr),]

cbPalette_viridis <- viridis(9,option = "D")
fit.ggplot_exp=data.frame(y=EXP_Log_fn(xdata$dist,coef(Expo_fit_result)[[1]],coef(Expo_fit_result)[[2]]),x=xdata$dist)
fit.ggplot_lomax_1=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[1]])[[2]],coef(all_t_lomax_est[[1]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[1]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_2=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[2]])[[2]],coef(all_t_lomax_est[[2]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[2]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_3=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[3]])[[2]],coef(all_t_lomax_est[[3]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[3]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_4=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[4]])[[2]],coef(all_t_lomax_est[[4]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[4]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_5=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[5]])[[2]],coef(all_t_lomax_est[[5]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[5]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_6=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[6]])[[2]],coef(all_t_lomax_est[[6]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[6]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_7=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[7]])[[2]],coef(all_t_lomax_est[[7]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[7]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_8=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[8]])[[2]],coef(all_t_lomax_est[[8]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[8]/4)^2)))),x=xdata$dist)
fit.ggplot_lomax_9=data.frame(y=LOMAX_Log_fn(xdata$dist,coef(all_t_lomax_est[[9]])[[2]],coef(all_t_lomax_est[[9]])[[1]],1/(t_expo_est*(t_expo_est/((t_d[9]/4)^2)))),x=xdata$dist)

ggplot(data=xdata,aes(x=dist,y=log(wcorr)))+
  geom_point()+
  xlim(0,1)+
  ylim(-15,-5)+
  geom_line(data=fit.ggplot_exp,aes(x=x,y=y),color=cbPalette_viridis[1])+
  geom_line(data=fit.ggplot_lomax_1,aes(x=x,y=y),color=cbPalette_viridis[1])+
  geom_line(data=fit.ggplot_lomax_2,aes(x=x,y=y),color=cbPalette_viridis[2])+
  geom_line(data=fit.ggplot_lomax_3,aes(x=x,y=y),color=cbPalette_viridis[3])+
  geom_line(data=fit.ggplot_lomax_4,aes(x=x,y=y),color=cbPalette_viridis[4])+
  geom_line(data=fit.ggplot_lomax_5,aes(x=x,y=y),color=cbPalette_viridis[5])+  
  geom_line(data=fit.ggplot_lomax_6,aes(x=x,y=y),color=cbPalette_viridis[6])+
  geom_line(data=fit.ggplot_lomax_7,aes(x=x,y=y),color=cbPalette_viridis[7])+
  geom_line(data=fit.ggplot_lomax_8,aes(x=x,y=y),color=cbPalette_viridis[8])+
  geom_line(data=fit.ggplot_lomax_9,aes(x=x,y=y),color=cbPalette_viridis[9])+
  labs(y = "log weighted LD")+
  labs(x = "Genetic Distance in cM")

ggplot(data=Mean_time_est, aes(x=Model, y=Mean_time,fill=Model)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c(cbPalette_viridis[1],cbPalette_viridis))


