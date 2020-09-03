FUN <-function(A,doses){
  b <- A[1] + A[2]*exp((doses-A[5])^2*(-A[6]))*(1/(1+exp(-(doses-A[3])/A[4])))
  return(b)
}

dFUN <- function(A,d){
  a = -2*A[6]*(d-A[5]) + (1/A[4])*exp(-(1/A[4])*(d-A[3]))/(1  + exp(-(1/A[4])*(d-A[3])))
  return (a)
}


d2FUN <- function(A,d){
  a = -2*A[6] - (1/A[4])^2*exp(-(1/A[4])*(d-A[3]))/((1  + exp(-(1/A[4])*(d-A[3]))))^2
  return (a)
}

library(ToxicR)

set.seed(893223)

D <-c(rep(seq(0,1.0,1/4),each=5))
mean <- 2.3  + 10/(1+exp(-(D-0.77)*5))*(1/(1+exp(-(0.75-D)*3)))

Y <- mean + rnorm(length(mean),0,0.2)
Q <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="FUNL",distribution = "normal",fit_type = "laplace")
R <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="FUNL",distribution = "normal",fit_type = "mcmc")
S <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="hill",distribution = "normal",fit_type = "mcmc")
Z <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="exp-5",distribution = "normal",fit_type = "laplace")
U <- single_continuous_fit(as.matrix(D),as.matrix(Y),sstat = F,BMR = 1.0 ,model_type="exp-5",distribution = "normal",fit_type = "mcmc")
#system.time({B <-test_skn(as.matrix(Y), as.matrix(D))})
B <- Q$fitted_model$parameters

plot(D,Y,pch=16)

doses <- seq(0,1,0.01)

lines(doses,FUN(B,doses),col=2,lwd=3)

