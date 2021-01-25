# Purpose: Test data visualization part for ToxicR package
# Author: Sooyeong Lim
# Last Updated Date: 01/24/2021


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

# Case 1: Dichotomous - MCMC Fitting
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "mcmc")
# Test 1. Dichotomous MA Clevland Plot
.cleveland_plot.BMDdichotomous_MA(A)
# Test 2. Dichotomous MA Density Plot 
.plot.density.BMDdichotomous_MA_MCMC(A)


# Case 2: Dichotomous - laplace
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "laplace")
# Test 1. Dichotomous MA Clevland Plot
.cleveland_plot.BMDdichotomous_MA(A)
# Test 2. Dichotomous MA Density Plot 
#-- This part should be double checked. Laplace fitting is not consistent with other fitting model
.plot.density.BMDdichotomous_MA_maximized(A)




# Sample Data - Continous Example

D <-c(rep(seq(0,1.0,1/4),each=4))
mean <- 2.3  + 10/(1+exp(-(D-0.60)*8))*(1/(1+exp(-(0.99-D)*13)))
Y <- mean + rnorm(length(mean),0,0.7)


# Case 1: Dichotomous - MCMC Fitting

# Question- continous case why there are 13 models?
A<-ma_continuous_fit(D,Y,fit_type="mcmc",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd',
                    model_list = c("hill", "exp-3", "exp-5", "power", "FUNL"),
                     distribution_list = c(rep("normal",5)))

# Test 1. Dichotomous MA Clevland Plot
.cleveland_plot.BMDcontinous_MA(A)
.plot.BMDcontinuous_MA(A)
# Bit weird result FUNL- Almost 1 other cases is 0;
.plot.density.BMDcontinous_MA_MCMC(A)


