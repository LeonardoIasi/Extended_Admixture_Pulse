---
title: "Fitting Simple and extended admixture pulses to segment data"
author: "Leonardo N. M Iasi"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Fitting the simple and extended pulse to segment length data

In this documentation we will discuss two ways of modeling the length distribution of introgressed segments after a certain time t. The simple pulse model assumes that all the introgressed materia entered in just one generation (i.e. the duration of gene flow was one gneneration) a time $t_m$ ago. The extended pulse relaxes this assumption. here, the introgressed material can enter over an extended period of time with a mean time of that period being $t_m$. The duration ($t_d$) is defined by the parameter k, wherease $t_d = 4 t_m k^{(-\frac{1}{2})}$. 

We will follow an example of an extended period of gene flow between Neandertal and Non-Africans around 50,000 years ago. We use simulated data generated using msprime coalescent simulations. The simulation used an emperical map to simulate realistic recombination patterns. In this simulation $t_m = 1,500$ generations ago and $t_d = 2000$ generations. The simulated segments are given in the Example_seg.txt file.


# Fitting the simple to segment length data

First, we want to fit the simple pulse model, which assumes implicitly a duration $t_d = 1$. The simple pulse will give us the estimate for the mean time of admixture $t_m$.

## Parameters

This are the parameters used for the function which fitts the simple pulse model.

### The path to the segment file


The segment file is read in using the [input] parameter (file_path). The file must at least contain one column giving the length of each unique segment in cM. The column must be named [length_cM].

```{r}
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

input <- "Example_seg.txt"
```


### truncation

Somtimes it is necessery to truncate the segment length distribution, especially if the segments are inferred you might want to exclude extremly short segments (e.g. smaller than 0.05 cM) since they are usually inferred with low confidence. You might also exclude very long segments wich might be wrongly inferred by some approaches where two segments close by are accedently joined to one long segment. The upper length cutoff is highly dependent on your scenario you are interrested in. For Neandertal gene flow we do not expect segments to be longer then 1 cM.

Parameter (truncation = TRUE/FALSE) to define if you want to exclude segments given a ceartain threshold lengtn in cM from the fitting (defined by upper and lower trunc parameters)

```{r}
truncation <-  T
```

### upper and lower thresholds

Only used if truncation = T. Numeric parameter to define lower (lower_trunc) and upper (upper_trunc) cutoff for the segment length given in cM 

```{r}
lower_trunc <- 0.05	
upper_trunc <- 1
```

### tm upper and lower boundries

To limit the searchspace, especially for the extended pulse model fit you must give an upper and lower boundry for the $t_m$ parameter. In our example we know from Archaeological records that no Neandertal remains were found to be younger than 35,000 years so we can set a very loose lower boundry at 100 generations ago (3000 years ago) and an upper boundry 5000 generations ago, which is a roughly double the estimated  split time between Africans and Non-Africans. You need to adjust these boundrys to your model organism. At least give a lower boundry at 1 generation and a finite value for the upper since the time of gene flow can not be negative.

```{r}
tm_lower <-  100
tm_upper <-  5000
```

## Running the function

Now we defind all our parameters. All we need to doo now is executing the function.

```{r}
fit_simple_pulse=function(input,truncation=F,lower_trunc=NA,upper_trunc=NA,tm_lower,tm_upper){
  Segments <- read.table(Segments,header = T)
  
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
    return(list(res,data.frame(t_m=t_m, ll_exp = -res$value,lower_trunc=lower_trunc,upper_trunc=upper_trunc),l))
  }, silent=F)
  return(list(NA,data.frame(t_m=NA, ll_exp=NA ,lower_trunc=lower_trunc,upper_trunc=upper_trunc),l))
}

SP_Fit <- fit_simple_pulse(input = input,truncation = truncation,lower_trunc =lower_trunc,upper_trunc = upper_trunc,
                           tm_lower = tm_lower,tm_upper = tm_upper)

tm_SP_est <- SP_Fit[[2]]$t_m
ll_SP_est <- SP_Fit[[2]]$ll_exp
lower_trunc_tm_SP_est <- SP_Fit[[2]]$lower_trunc
upper_trunc_tm_SP_est <- SP_Fit[[2]]$upper_trunc

```

## Visualize the fit to the data


# Fitting the extended to segment length data
 
We use the same setup and data as before but now we want to fit the extended pulse which has one additional parameter $k$. This parameter is related to the duration of admixture such that $t_d = 4 t_m k^{(-\frac{1}{2})}$. 

## Parameters

The extenmded pulse function takes the same parameter with two additional ones setting the boundries for k.

### k upper and lower boundries

This is only needed for the extended pulse model. The parameter k is defined for a parameter space between 2 (continuous gene flow) and infinity (pulse gene flow). Realistically the parameter will not go to infinity and will fit in a one generation pulse with a high number at least  bigger than 100. here we choose 1e8. Our lower value is 2, since k is not defined for a value smaller than 2.
```{r}
k_lower <- 2
k_upper <- 1e8
```


```{r, echo=FALSE, results='asis'}
knitr::kable(head(mtcars, 10))
```

