library(ToxicR)

set.seed(893223)

D <-c(rep(seq(0,1.0,1/6),each=4))
#D <- runif(100)
mean <- 2.3  + 10/(1+exp(-(D-0.60)*8))*(1/(1+exp(-(0.99-D)*13)))

Y <- mean + rnorm(length(mean),0,3.15)
D = D*10
a <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMD_TYPE = "sd" ,BMR = 1 ,model_type="FUNL",distribution = "normal-ncv",fit_type = "laplace")
b <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMD_TYPE = "sd" ,BMR = 2 ,model_type="hill",distribution = "normal",fit_type = "laplace")
c <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMD_TYPE = "sd" ,BMR = 2 ,model_type="exp-5",distribution = "normal",fit_type = "laplace")
d <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMD_TYPE = "sd" ,BMR = 2 ,model_type="exp-3",distribution = "normal",fit_type = "laplace")
e <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMD_TYPE = "sd" ,BMR = 2 ,model_type="power",distribution = "normal",fit_type = "laplace")

fit<-ma_continuous_fit(D,Y,fit_type="mcmc",   samples=25000,burnin=2500,BMR=1.0,BMD_TYPE = "sd")
fit<-ma_continuous_fit(D,Y,fit_type="laplace",samples=25000,burnin=2500,BMR=1.0,BMD_TYPE = "sd")


#system.time({fit2<-ma_continuous_fit(D,Y,fit_type="laplace",BMR=2.0)})

#plot(fit2)
plot(fit)

D <-c(rep(seq(0,1.0,1/4),each=4))
mean <- 2.3  + 10/(1+exp(-(D-0.60)*5))*(1/(1+exp(-(0.99-D)*10)))

Y <- mean + rnorm(length(mean),0,1.5)
#Q <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="FUNL",distribution = "normal",fit_type = "laplace")

fit<-ma_continuous_fit(D,Y,fit_type="mcmc",samples=25000,burnin=2500,BMR=2.0)
fit2<-ma_continuous_fit(D,Y,fit_type="laplace",BMR=2.0)

plot(fit2)
plot(fit)

plot(D,Y,pch=16)

doses <- seq(0,1,0.01)

lines(doses,FUN(B,doses),col=2,lwd=3)

