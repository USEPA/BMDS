library(ToxicR)

mData <- matrix(c(0, 1,10,
                  0.3, 4,10,
                  1, 4,10,
                  4, 7,10),nrow=4,ncol=3,byrow=T)


mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)

prior <- create_prior_list(normprior(	0,	10,	-18,	18),
                           lnormprior(	0.693147180559945,	0.01,	0.2,	20),
                           lnormprior(	0,	2,	0,	1e4))

G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "laplace")
system.time({B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mle")})

system.time({A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3], fit_type = "mle")})


G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")
G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mle")


G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma", fit_type = "laplace")
G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma", fit_type = "mle")


A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "laplace")

D = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "laplace")
D = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "mle")

E = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "laplace")
E = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "mle")

H = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "laplace")
H = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "mle")

I = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "laplace")
I = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "mle")

J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "qlinear",fit_type = "laplace")
J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "qlinear",fit_type = "mle")

J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit",fit_type = "laplace")
J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit",fit_type = "mle")

J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "laplace")
J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "mle")



B = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "laplace")
plot(G)

plot(K[1,],ylim=c(0,1))
for (ii in 1:2000){
  lines(K[ii,])
}
plot(B$Individual_Model_5)


mData <- matrix(c(0, 0,5,
                  0.1, 2,5,
                  0.3, 3,5,
                  1, 4,5,
                  10,5,5),nrow=5,ncol=3,byrow=T)
                  
mData <- matrix(c(0, 0,50,
                  1, 1,50,
                  3, 2,50,
                  7, 5,50,
                  10, 20,50),nrow=5,ncol=3,byrow=T)

mData <- matrix(c(0, 0,10,
                  2, 0,10,
                  4, 0,10,
                  10, 10,10),nrow=4,ncol=3,byrow=T)

mData <- matrix(c(0, 0,50,
                  62.5, 9,50,
                  125, 8,50,
                  250, 14,50),nrow=4,ncol=3,byrow=T)


A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "mle")
plot(B,pch=16,ylab="Probability: adenomas + carcinomas",xlab="1-Bromopropane (ppm)")
plot(C,pch=16,ylab="Probability: adenomas + carcinomas",xlab="1-Bromopropane (ppm)")

cdf_rows = nrow(C$cdf)

cdf_rows = nrow(fit$cdf)s
a = fit$cdf[,1]
b = fit$cdf[,2]
b = b[(!is.infinite(a))*(!is.na(a))*(!is.nan(a)) == 1]
a = a[(!is.infinite(a))*(!is.na(a))*(!is.nan(a)) == 1]
ta = a
tb = ta*0
a = c(min(a)*0.9,a,max(a)*1.1)
b = c(0,b,1.0)
cdf_est <- splinefun(a,b)
h = 1e-6
for(i in 1:length(ta)){
  tb[i] = (cdf_est(ta[i]+h)-cdf_est(ta[i]))/h
}
  
tb = tb/(max(tb))*0.4*0.4
dens
polygon(ta,tb,col = alphablend(col="red",0.2),border =alphablend(col="red",0.2))


lines(doses,lq)
