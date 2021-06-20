# Test examples for every plot functions
# Author: Sooyeong Lim
# Date: 06/20/21
library(ToxicR)
library(ggridges)
library(forcats)
library(ggplot2)
library(dplyr)
library(readr)




# Example 1- Dichotomous Example

mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)

# Single model dichotomous fitting

A_single_mcmc<-single_dichotomous_fit(mData[,1],mData[,2],mData[,3], model_type="hill",fit_type="mcmc")
# Updates
# 1. y axis is scaled based on data input
# 2. Color of density plot is updated to match with continous case
plot(A_single_mcmc)
.plot.BMDdich_fit_MCMC(A_single_mcmc)


A_single_laplace = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type="hill",fit_type = "laplace")
# Update
# 1. y axis is scaled based on data input
plot(A_single_laplace)
.plot.BMDdich_fit_maximized(A_single_laplace)


# Model Average dichotomous fitting
# Dichotomous - Model Average
# Case 1: Dichotomous - MCMC Fitting
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "mcmc")

# Base plot
plot(A)
# Base plot - MA density curve seems little bit odd
# 1. Alpha level adjust based on the posterior probabilities
# 2. Density plot color matching 
# 3. Re-scale y axis based on the data
.plot.BMDdichotomous_MA(A)
# Test 1. Dichotomous MA Cleveland Plot

# Hold- Update cleveland plot
# 1. Posterior probability check 
.cleveland_plot.BMDdichotomous_MA(A)
# Test 2. Dichotomous MA Density Plot
# 1. Color adjust based on John's comment - Red for MA and Grey for other models
.plot.density.BMDdichotomous_MA_MCMC(A)



# Case 2: Dichotomous - laplace MA
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "laplace")
# Test 1. Dichotomous MA Clevland Plot-- For laplace case/works
.cleveland_plot.BMDdichotomous_MA(A)

# This part needs to be double checked & density part should be re-derived
# - Code / Need to assign dynamically based on how many models are there
.plot.BMDdichotomous_MA(A)






# Example 2- Continuous case example

D <-c(rep(seq(0,1.0,1/4),each=4))
mean <- 2.3  + 10/(1+exp(-(D-0.60)*8))*(1/(1+exp(-(0.99-D)*13)))
Y <- mean + rnorm(length(mean),0,0.7)


# Continuous example with sufficient statistic==F 

A<- single_continuous_fit(as.matrix(D),as.matrix(Y),model_type = "exp-5",
                          distribution="normal-ncv",fit_type = "mcmc",sstat = F)

plot(A)
.plot.BMDcont_fit_MCMC(A)

A<- single_continuous_fit(as.matrix(D),as.matrix(Y),model_type = "exp-5",
                          distribution="normal-ncv",fit_type = "laplace",sstat = F)

plot(A)
.plot.BMDcont_fit_maximized(A)


# Fitting example should be much simpler - 
A<-ma_continuous_fit(D,Y,fit_type="mcmc",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd')


# This part's alpha part needs to be fixed as of dichotomous case;

# Updates
# 1. Density kernel's n=512
# 2. Alpha level is adjusted based on the model's posterior probability
plot(A)
.plot.BMDcontinuous_MA(A)


# Update based on prior probability for continous case 
.cleveland_plot.BMDcontinous_MA(A)
.plot.density.BMDcontinous_MA_MCMC(A)
