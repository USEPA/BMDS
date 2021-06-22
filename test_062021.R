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

# - Code / Need to assign dynamically based on how many models are there
plot(A)
.plot.BMDdichotomous_MA(A)






# Example 2- Continuous case example

D <-c(rep(seq(0,1.0,1/4),each=4))
mean <- 2.3  + 10/(1+exp(-(D-0.60)*8))*(1/(1+exp(-(0.99-D)*13)))
Y <- mean + rnorm(length(mean),0,0.7)


# Continuous example with sufficient statistic==F 

A<- single_continuous_fit(as.matrix(D),as.matrix(Y),model_type = "exp-5",
                          distribution="normal-ncv",fit_type = "mcmc",sstat = F)

plot(A)
# Need to check this part - 06/21/21 4:30PM 

.plot.BMDcont_fit_MCMC(A)

A<- single_continuous_fit(as.matrix(D),as.matrix(Y),model_type = "exp-5",
                          distribution="normal-ncv",fit_type = "laplace",sstat = F)

plot(A)
.plot.BMDcont_fit_maximized(A)




model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))



# Fitting example should be much simpler - 
A<-ma_continuous_fit(D,Y,fit_type="mcmc",samples=25000,burnin=2500,BMR=0.1,BMD_TYPE='sd', model_list=model_list)


# This part's alpha part needs to be fixed as of dichotomous case;

# Updates
# 1. Density kernel's n=512
# 2. Alpha level is adjusted based on the model's posterior probability
plot(A)
.plot.BMDcontinuous_MA(A)

# Update based on prior probability for continuous case 
.cleveland_plot.BMDcontinous_MA(A)
.plot.density.BMDcontinous_MA_MCMC(A)







# Sufficient statistics test


M2           <- matrix(0,nrow=5,ncol=4)
colnames(M2) <- c("Dose","Resp","N","StDev")
M2[, 1]      <- c(0,50, 100, 150, 200)
M2[, 2]      <- c(10, 20 , 30, 40 ,50)
M2[, 3]      <- c(100, 100, 100, 100, 100)
M2[, 4]      <- c(3, 4, 5, 6, 7)


# I think sstat=F or T doesn't matter since it reads # of columns
c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],sstat=T,BMD_TYPE="sd",BMR=1, 
                          distribution = "normal",fit_type="mle",model_type = "power")

# Adjust size of the interval bar;


plot(c)
#Color and style needs to be matched
.plot.BMDcont_fit_maximized(c,qprob=0.05)


c_mcmc = single_continuous_fit(M2[,1,drop=F],M2[,2:4],sstat=T,BMD_TYPE="sd",BMR=1, 
                          distribution = "normal",fit_type="mcmc",model_type = "power")

plot(c_mcmc)
#Color and style needs to be matched
.plot.BMDcont_fit_MCMC(c_mcmc,qprob=0.05)







doses	<- c(0,	0,	0,	0,	0.156,	0.156,	0.156,	0.3125,	0.3125,	0.3125,
           0.625,	0.625,	0.625,	1.25,	1.25,	1.25,	2.5,	2.5,	2.5,	5,5,
           5,	5,	10,	10,	10,	10,	20,	20,	20,	20) 

PFOA_Liver <- read_table2("PFOA_Liver.txt", 
                          col_names = FALSE)


doses	<- c(0,	0,	0,	0,	0.156,	0.156,	0.156,	0.3125,	0.3125,	0.3125,
           0.625,	0.625,	0.625,	1.25,	1.25,	1.25,	2.5,	2.5,	2.5,	5,5,
           5,	5,	10,	10,	10,	10,	20,	20,	20,	20) 


# Model average fitting 

model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))



# 05/28 SL Try to focus on the issue here
# 1. CI Band seems odd
# 2. BMD Density plot looks odd 
# 3. Color Theme update for giving user more option

# I think it is useful that which model is dominant in terms of the line- color
# Range max-min part should be updated... 
# BMDU is out of bound


# How about the log scale?

BB <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "mcmc",BMR = 1,model_list = model_list)

# X - axis is bit odd for log10() case 
# For model average why it shows only one model?

plot(BB)



# Question not sure why BMD is superimposed? 
.plot.BMDcontinuous_MA(BB, qprob=0.05)

.plot.BMDcontinuous_MA(BB, qprob=0.05)+scale_x_log10()

# Update based on prior probability for continous case 
.cleveland_plot.BMDcontinous_MA(BB)
.plot.density.BMDcontinous_MA_MCMC(BB)




