#set.seed(12345)
#library(BMDS) # Uncomment this line if BMDS is not already manually loaded in R console
#library(Rcpp)

log_likelihood <- function(Y,mu,var){
  M_PI = pi
  returnV = 	Y[,3]*log(1/sqrt(2.0*M_PI)) - (Y[,3] / 2.0)*log(var) - 
              (1 / (2.0*var))*((Y[,3] - 1)*Y[,3]^2 +
              Y[,3]*(Y[,1] - mu)^2);
  return(sum(returnV))
}


############################################################################
####################     Data set: Continuous3.dax     #####################
############################################################################
library(ToxicR)
c3ss <- T # Summarized dose-response data

c3Dat <- matrix(0,nrow=4,ncol=4)
colnames(c3Dat) <- c("Dose","Resp","N","StDev")
c3Dat[, 1] <- c(0,10,28,33)*1000
c3Dat[, 2] <- c(6.1,8.1,11.5,10.6) 
c3Dat[, 3] <- c(5,5,5,5)
c3Dat[, 4] <- c(0.93,1.94,5.73,2.61)

system.time({ BB <- single_continuous_fit(as.matrix(c3Dat[, 1]),c3Dat[, 2:4],model_type = "exp-5",degree = 3,BMD_TYPE = "sd",BMR= 1,
                    distribution="normal",fit_type = "laplace",isFast=T) })


AA <- ma_continuous_fit(as.matrix(c3Dat[, 1]),c3Dat[, 2:4],BMD_TYPE = "rel",BMR= 0.1,
                        fit_type = "mcmc")

#print("Running EXP-3 Opt1")

BB <- single_continuous_fit(as.matrix(c3Dat[, 1]),c3Dat[, 2:4],model_type = "polynomial",BMD_TYPE = "sd",BMR= 1,
                            distribution="normal",fit_type = "laplace")

#BB$covariance
#BB$parameters

AA <- ma_continuous_fit(as.matrix(c3Dat[, 1]),c3Dat[, 2:4],BMD_TYPE = "rel",BMR= 0.1,
                            fit_type = "laplace")
plot(AA)

k <- seq(0,16.6,by=0.05)
values <- matrix(BB$parameters[1:3],ncol=3,nrow=length(k),byrow=T)
values[,3] = k
likelihood <- 0*k
dim(values)
for(ii in 1:length(k)){
  mean = cont_power_f(c(values[ii,]),c3Dat[, 1])
  likelihood[ii] = log_likelihood(c3Dat[, 2:4],mean,1)
}

