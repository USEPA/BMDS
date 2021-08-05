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

prior <- create_prior_list(normprior(0,10,-100,100),
                           lnormprior(0,1,0,100),
                           lnormprior(0,1,0,100));

a <- create_dichotomous_prior(prior,"gamma")
b <- create_dichotomous_prior(prior,"weibull")
prior2 <- create_prior_list(normprior(0,10,-100,100),
                            normprior(0,1 ,-100,100),
                            lnormprior(0,1,0,100));

c <- create_dichotomous_prior(prior2,"log-logistic")

vi <- list(a,b,c)

library(ToxicR)

mData <- matrix(c(0,	8,	50,
                  209.8,	17,	50,
                  444.6,	26,	50,
                  978.1,	42,	50),nrow=4,ncol=3,byrow=T)
#Q  = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_list = vi, fit_type = "mcmc")
#Q1 = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3], fit_type = "mcmc")
#plot(Q) + scale_x_continuous(trans="sqrt")


system.time({R = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",degree = 3, fit_type = "mcmc")})
S = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",degree = 3, fit_type = "laplace")

R = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],prior = q,degree = 2,fit_type = "laplace")
S = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",degree = 2,fit_type = "laplace")

plot(R)

R = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",degree = 2,fit_type = "mcmc")
