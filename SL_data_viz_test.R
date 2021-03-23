# Purpose: Test data visualization part for ToxicR package
# Author: Sooyeong Lim
# Last Updated Date: 03/03/2021

# Load packages
library(ToxicR)
library(ggridges)
library(forcats)
library(ggplot2)
library(dplyr)


# Sample Data - Dichotomous Example
mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)


# Single model fitting- MCMC


# too much white space- this needs to be adjusted automatically

A_single_mcmc<-single_dichotomous_fit(mData[,1],mData[,2],mData[,3], model_type="hill",fit_type="mcmc")
.plot.BMDdich_fit_MCMC(A_single_mcmc)

# Single model fitting- Laplace
A_single_laplace = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type="hill",fit_type = "laplace")
.plot.BMDdich_fit_maximized(A_single_laplace)
# Base density plot is required? 
 


# Dichotomous - Model Average
# Case 1: Dichotomous - MCMC Fitting
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "mcmc")
# Test 1. Dichotomous MA Clevland Plot
.cleveland_plot.BMDdichotomous_MA(A)
# Test 2. Dichotomous MA Density Plot - Update for base-color later
.plot.density.BMDdichotomous_MA_MCMC(A)
### - Need to apply John's idea 

#Base plot - MA density seems little bit odd
.plot.BMDdichotomous_MA(A)


# Case 2: Dichotomous - laplace
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "laplace")
# Test 1. Dichotomous MA Clevland Plot-- For laplace case/works
.cleveland_plot.BMDdichotomous_MA(A)

# This part needs to be double checked & density part should be re-derived
# - Code / Need to assign dynamically
.plot.BMDdichotomous_MA(A)






# Question to Matt.
# This part should be double checked. Laplace fitting is not consistent with other fitting model

# This one shows the CDF of BMD 
A$Fitted_Model_1$bmd_dist
# Base plot needs to be fixed too. It has issue with Model- MA representation.
plot(A)
# This is point estimates 50% , 5%, 95%...
# How can we set MA sample desnity plot from the output data from here? 
A$Fitted_Model_1$bmd



# Sample Data - Continous Example

D <-c(rep(seq(0,1.0,1/4),each=4))
mean <- 2.3  + 10/(1+exp(-(D-0.60)*8))*(1/(1+exp(-(0.99-D)*13)))
Y <- mean + rnorm(length(mean),0,0.7)
plot(D,Y)


# Matt's code update
#plot(doses,sim_data[1,])


# Case 1: Dichotomous - MCMC Fitting
# Question to Maatt continous case why there are 13 models?
# Different Prior's- Should we show them separately?

# Continous single fit is required to be tested
k<-single_continuous_fit(D, Y, model_type = "hill", fit_type = "mcmc", prior = "default", 
                      BMD_TYPE = "sd", sstat = T, BMR = 0.1, point_p = 0.01, distribution = "normal-ncv", 
                      alpha = 0.05, samples = 21000, degree = 2, burnin = 1000)

View(single_continuous_fit)


# What I need to change is .. 1. Too much white space in the figure/ plot
# Fitting example should be much simpler - 
A<-ma_continuous_fit(D,Y,fit_type="mcmc",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd')

plot(A)


# There will be 13 models in this case- Let's fix something easier

# Condition might be added - if it is less than ... 0.05


A$posterior_probs

# Test 1. Dichotomous MA Clevland Plot
.cleveland_plot.BMDcontinous_MA(A)
# Bit weird result FUNL- Almost 1 other cases is 0;
.plot.density.BMDcontinous_MA_MCMC(A)
# Base object should be updated by using Shiny App
A$posterior_probs
#It is dominated by the FUNL model. While the other 4 models are minor 

# This continous baseplot should be updated
.plot.BMDcontinuous_MA(A)





A<-ma_continuous_fit(D,Y,fit_type="laplace",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd',
                     model_list = c("hill", "exp-3", "exp-5", "power", "FUNL"),
                     distribution_list = c(rep("normal",5)))
plot(A)
