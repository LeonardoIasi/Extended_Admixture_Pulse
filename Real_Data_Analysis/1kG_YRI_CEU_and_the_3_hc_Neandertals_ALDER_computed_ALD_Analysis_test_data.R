wcorr= exp(-(seq(0.05,3,0.001)/100)/(1/1500))
xdata_test=data.frame(dist=seq(0.05,3,0.001),wcorr=wcorr)
plot(xdata_test$dist,xdata_test$wcorr)


Expo_fit_result_test <- Fit_Exp_fn(Data = xdata_test,affine = T)
t_expo_est=1/coef(Expo_fit_result_test)[[2]]

t_expo_est=t_expo_est
t_d=c(1,100,200,400,800,1000,1500,2000,2500)
all_lomax_fit <- list()
all_t_lomax_est <- c()
for(t in 1:length(t_d)){
  fixed_w= 1/(t_expo_est*(t_expo_est/((t_d[t]/4)^2)))
  
  Lomax_fit_result <- Fit_Lomax_fn(Data = xdata_test,fixed_w = fixed_w)
  coef(Lomax_fit_result)
  t_lomax_est=1/coef(Lomax_fit_result)[[1]]
  all_t_lomax_est <- c(all_t_lomax_est,t_lomax_est)
  all_lomax_fit[[t]] <- Lomax_fit_result
}



AIC_result <- AIC(Expo_fit_result_test,all_lomax_fit[[1]],all_lomax_fit[[2]],all_lomax_fit[[3]],all_lomax_fit[[4]],
                  all_lomax_fit[[5]],all_lomax_fit[[6]],all_lomax_fit[[7]],all_lomax_fit[[8]],all_lomax_fit[[9]])


Model_RSS <- c()
for(i in list(Expo_fit_result_test,all_lomax_fit[[1]],all_lomax_fit[[2]],all_lomax_fit[[3]],all_lomax_fit[[4]],
               all_lomax_fit[[5]],all_lomax_fit[[6]],all_lomax_fit[[7]],all_lomax_fit[[8]],all_lomax_fit[[9]])){
  RSS <- sum(residuals(i)^2) 
  Model_RSS <- c(Model_RSS,RSS)
}

Mean_time_est <- data.frame(Model=c("Exp","Lomax(td=1)","Lomax(td=100)","Lomax(td=200)","Lomax(td=400)"
                                    ,"Lomax(td=800)","Lomax(td=1000)","Lomax(td=1500)","Lomax(td=2000)","Lomax(td=2500)"),
                            Mean_time=c(1/coef(Expo_fit_result_test)[[2]],1/coef(all_lomax_fit[[1]])[[1]],1/coef(all_lomax_fit[[2]])[[1]]
                                        ,1/coef(all_lomax_fit[[3]])[[1]],1/coef(all_lomax_fit[[4]])[[1]],1/coef(all_lomax_fit[[5]])[[1]]
                                        ,1/coef(all_lomax_fit[[6]])[[1]],1/coef(all_lomax_fit[[7]])[[1]],1/coef(all_lomax_fit[[8]])[[1]]
                                        ,1/coef(all_lomax_fit[[9]])[[1]]))


pdf("Progress_figures_test_data.pdf")
plot(xdata_test$dist,(xdata_test$wcorr),pch=19,xlim=c(0,0.5))
lines(xdata_test$dist,(predict(Expo_fit_result_test)),col="red", lwd=2)
for(t in 1:length(t_d)){
  lines(xdata_test$dist,(predict(all_lomax_fit[[t]])),col="green", lwd=2)
}
barplot(Mean_time_est$Mean_time,names.arg = Mean_time_est$Model,las=2,main="mean time est")
barplot(Model_RSS,names.arg = Mean_time_est$Model,las=2,main="RSS")
barplot(AIC_result$AIC+44000,names.arg = Mean_time_est$Model,las=2,main="AIC")
dev.off()