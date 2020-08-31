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

D <- rep(seq(0,1.0,1/4),each=5)
mean <- 2.3  + 2/(1+exp(-(D-0.37)*5))*(1/(1+exp(-(0.75-D)*3)))

Y <- mean + rnorm(length(mean),0,0.5)
system.time({B <-test_skn(as.matrix(Y), as.matrix(D))})


plot(D,Y,pch=16)

doses <- seq(0,1,0.01)


FUN(a,0.08636607) - FUN(a,0)
FUN(b,0.08811084) - FUN(b,0)

d = 10
for (ii in 1:7){
  d = d - dFUN(B,d)/d2FUN(B,d)
}

lines(doses,FUN(B,doses),col=2,lwd=3)

