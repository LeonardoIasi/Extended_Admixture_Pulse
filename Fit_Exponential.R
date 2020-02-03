suppressPackageStartupMessages({
  library("DEoptim")
  library("MASS")
})
#!/usr/bin/Rscript
args <- commandArgs(TRUE)

######### Fit exponential and/or lomax distribution to ROLLOFF output ############
# rexpfit.r is a program that can be used to fit an exponential distribution to the output of ROLLOFF. This program uses the nls function in R to determine the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model. It also uses the DEoptim package to pick the initial values of all the parameters (Katharine Mullen, David Ardia, David Gil, Donald Windover, James Cline (2011). DEoptim: An R Package for Global Optimization by Differential Evolution. Journal of Statistical Software, 40(6), 1-26. URL http://www.jstatsoft.org/v40/i06/.).
#More information on these functions and packages can be found at-
#http://stat.ethz.ch/R-manual/R-patched/library/stats/html/nls.html
#http://cran.r-project.org/web/packages/DEoptim/index.html

#####################################################################


input <- as.character(args[1]) 		# output from rollof
Snake_output_One <- as.character(args[2])  		# output of expfit log file
Snake_output_Two <- as.character(args[3]) # output of expfit plot of the decay
lval <- as.numeric(args[4])		# lower value of dist to use
hval <- as.numeric(1)		# higher value of dist to use  
col <- as.numeric(2)   # data column containing weighted correlation
plot <- as.logical(FALSE)           # plot the output = TRUE/FALSE 


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

### Fitting a Exopnential function ###

# Pick optimized parameters for model 1: A*exp(m1 * -x) + C
fm1_exp <- function(x) x[1]*exp(-x[2] * dist/100)+x[3]
fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,1,-1), upper=c(1, 5000,1), control=list(trace=FALSE))
par1_exp <- fm3_exp$optim$bestmem
# parameters for y ~ Ae-mt
names(par1_exp) <- c("A", "m","c")


# fit the model of Single exponential by triyng non-linear least square
Non_Linear_LS <- function(DEOptim_res){
  fit1_exp <- nls(wcorr ~ (A*exp(-m * dist/100)+c), start=DEOptim_res, control=list(maxiter=10000, warnOnly=TRUE))
  #estimated parameters
  A_exp_est <- coef(fit1_exp)[1]
  Lambda_est <- coef(fit1_exp)[2]  	# rate of decay of exponential
  C_exp_est <- coef(fit1_exp)[3]
  #res.conf <- cbind(org=coefficients(fit1_exp),confint(fit1_exp))
  #return(list(A_exp_est,Lambda_est,C_exp_est,"no_error",fit=as.list(fit1_exp),res.conf))
  return(list(A_exp_est,Lambda_est,C_exp_est,"no_error",fit=as.list(fit1_exp)))
}

alternativeFunction <- function(xx){return(c(xx,"nls_error",xx))}
options(warn = 2)
End_Fit=Try_Non_Linear_LS <-  try(Non_Linear_LS(par1_exp),silent = FALSE)
if (inherits(Try_Non_Linear_LS, "try-error")) End_Fit=alternativeFunction(par1_exp)

A_exp_est <- as.numeric(End_Fit[1])
Lambda_est <- as.numeric(End_Fit[2])  	# rate of decay of exponential
C_exp_est <- as.numeric(End_Fit[3])

Optim_Params_exp <- as.numeric(c(A_exp_est,Lambda_est,C_exp_est))
Get_Residual_exp <-function(x) (wcorr-fm1_exp(x))
Get_RSS_Sum_exp <-function(x) sum((wcorr-fm1_exp(x))^2)
Expo_Result= paste("Exponential function Fitted with","\n"," A =",A_exp_est=format(A_exp_est, digits=8, width=0, justify="right"),"Lambda =",Lambda_est=format(Lambda_est, digits=8, width=0, justify="right"),"c =",C_exp_est=format(C_exp_est, digits=8, width=0, justify="right"),"error=",End_Fit[4],sep=" ")


#print output
outlog <- paste(Snake_output_One)
cat("Summary of fit:\n", file=outlog, append = FALSE)

# Convert distance in Morgans
capture.output(summary(End_Fit[5]), file = outlog, append = TRUE)
cat(paste("A, m, c, RSS_Expo, nls status:","\t",A_exp_est,"\t",Lambda_est,"\t",C_exp_est,"\t",RSS_Expo=Get_RSS_Sum_exp(Optim_Params_exp),"\t",End_Fit[4],"\n",sep = ""), file=outlog, append=TRUE)

# plot results
if (plot) {
  outpdf <- paste(Snake_output_Two)
  pdf(outpdf)
  plot(dist,wcorr, col= "black", pch=16,xlab = "Genetic Distance (cM)", ylab = "Weighted Correlation", main =(Expo_Result),cex.main= 1,sub=  paste("RSS=",Get_RSS_Sum_exp(Optim_Params_exp),sep = " "))
  lines(dist, fm1_exp(c(A_exp_est,Lambda_est,C_exp_est)), type="l", col= "red" )
  legend("topright",cex=0.75, pch=c(16,NA), lty=c(NA,1), col=c("black", "red"), legend=c("Admixture LD", "Exponential fit"))
  garbage <- dev.off()
}