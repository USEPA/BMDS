library(ToxicR)

set.seed(893223)

D <-c(rep(seq(0,1.0,1/4),each=5))
mean <- 2.3  + 10/(1+exp(-(D-0.99)*5))*(1/(1+exp(-(0.99-D)*3)))

Y <- mean + rnorm(length(mean),0,0.1)
Q <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="FUNL",distribution = "normal",fit_type = "laplace")

system.time({fit<-ma_continuous_fit(D,Y,fit_type="mcmc")})
sprintf("%1.2f",fit$ma_results$posterior_probs)
R <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="power",distribution = "normal",fit_type = "mcmc",samples = 300000,burnin =  30000)
S <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="hill",distribution = "normal",fit_type = "mcmc")
Z <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="exp-5",distribution = "normal",fit_type = "laplace")
U <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="exp-5",distribution = "normal",fit_type = "mcmc")
#system.time({B <-test_skn(as.matrix(Y), as.matrix(D))})
B <- Q$fitted_model$parameters

plot(D,Y,pch=16)

doses <- seq(0,1,0.01)

lines(doses,FUN(B,doses),col=2,lwd=3)

