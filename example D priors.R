library(ToxicR)

temp_prior = create_prior_list(normprior(0,1,-18,18),
                               lnormprior(log(1.6),1,0,20),
                               lnormprior(0,1,0,20))

weibl  = create_dichotomous_prior(temp_prior,"weibull")
weibl

temp_prior = create_prior_list(normprior(0,1,-18,18),
                               lnormprior(0,1,0,20))

logist = create_dichotomous_prior(temp_prior,"logistic") 

logist


temp_prior = create_prior_list(normprior(0,1,-18,18),
                               normprior(0,1,0,20),
                               lnormprior(0,1,0,18))

loglogist = create_dichotomous_prior(temp_prior,"log-logistic")
loglogist

#put all of the priors into a list
models <- list(weibl,logist,loglogist)

mData <- matrix(c(0,    0,100,
                  50,   5,100,
                  100, 30,100,
                  150, 65,100,
                  200, 90,100),nrow=5,ncol=3,byrow=T)

fit <- ma_dichotomous_fit(as.matrix(mData[,1]),as.matrix(mData[,2]),as.matrix(mData[,3]),fit_type = "mcmc",
                          model_list = models)