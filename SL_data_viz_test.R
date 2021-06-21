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
# too much white space- this needs to be adjusted automatically - Not the for the first Priority, but needs to be updated once 
# Continuous part is updated

A_single_mcmc<-single_dichotomous_fit(mData[,1],mData[,2],mData[,3], model_type="hill",fit_type="mcmc")

#Need to color Match with Continuous plot 
plot(A_single_mcmc)
.plot.BMDdich_fit_MCMC(A_single_mcmc)

# Single model fitting- Laplace
A_single_laplace = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type="hill",fit_type = "laplace")
.plot.BMDdich_fit_maximized(A_single_laplace)
# Base density plot is required? 
 


# Dichotomous - Model Average
# Case 1: Dichotomous - MCMC Fitting
A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "mcmc")
# Test 1. Dichotomous MA Cleveland Plot
.cleveland_plot.BMDdichotomous_MA(A)
# Test 2. Dichotomous MA Density Plot - Update for base-color later
.plot.density.BMDdichotomous_MA_MCMC(A)
### - Need to apply John's idea 

# Density seems bit odd... need to adjust continuous case
plot(A)
#Base plot - MA density curve seems little bit odd
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


# Continous example  03/28/2021

A<- single_continuous_fit(as.matrix(D),as.matrix(Y),model_type = "exp-5",
                            distribution="normal-ncv",fit_type = "mcmc",sstat = F)


.plot.BMDcont_fit_MCMC(A)

A<- single_continuous_fit(as.matrix(D),as.matrix(Y),model_type = "exp-5",
                          distribution="normal-ncv",fit_type = "laplace",sstat = F)

.plot.BMDcont_fit_maximized(A)


# Fitting example should be much simpler - 
A<-ma_continuous_fit(D,Y,fit_type="mcmc",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd')


# This part's alpha part needs to be fixed as of dichotomous case;
.plot.BMDcontinuous_MA(A)
plot(A)

# Test 2. Dichotomous MA Cleveland Plot

.cleveland_plot.BMDcontinous_MA(A)
# The model's 
# Bit weird result FUNL- Almost 1 other cases is 0;
.plot.density.BMDcontinous_MA_MCMC(A)
# Base object should be updated by using Shiny App
A$posterior_probs
#It is dominated by the FUNL model. While the other 4 models are minor 

# This continuous base plot should be updated
.plot.BMDcontinuous_MA(A)







A<-ma_continuous_fit(D,Y,fit_type="laplace",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd')
.plot.BMDcontinuous_MA(A)



