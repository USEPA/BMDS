library(ToxicR)




mData <- matrix(c(0,    0,100,
                  50,   5,100,
                  100, 30,100,
                  150, 65,100,
                  200, 90,100),nrow=5,ncol=3,byrow=T)


mData <- matrix(c(0,    5,100,
                  50,   3,100,
                  100, 33,100,
                  150, 45,100,
                  200, 55,100),nrow=5,ncol=3,byrow=T)

library(ToxicR)
library(ggplot2)
library(dplyr)
mData <- matrix(c(0,   2,10,
                  15,   4,16,
                  18,   6,14,
                  21,   18,20),nrow=4,ncol=3,byrow=T)
Q = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "mcmc")
plot(Q) + scale_x_continuous(trans="sqrt")

R = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",degree = 2,fit_type = "mcmc")
plot(R)

R = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",degree = 2,fit_type = "mcmc")
