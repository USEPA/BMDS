# Test examples for every plot functions
# Author: Sooyeong Lim
# Date: 06/20/21

library(ToxicR)
library(ggridges)
library(forcats)
library(ggplot2)
library(dplyr)



# Example 1- Dichotomous Example

mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)

A_single_mcmc<-single_dichotomous_fit(mData[,1],mData[,2],mData[,3], model_type="hill",fit_type="mcmc")

#Need to color Match with Continuous plot 
plot(A_single_mcmc)
.plot.BMDdich_fit_MCMC(A_single_mcmc)
