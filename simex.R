library(actuar)
library(ToxicR)



x <- rinvgauss(1000, 10, shape = 200)

doses <- rep(c(0,25,50,75,100),each=5)
testd <- seq(0,100,1)

parms <- c(4,5.02,70,3.3)
parms <- c(4,5.02,40,1.3)
parms <- c(10.88,-3.02,25,3)

Exp5 <- rbind(c(481,0.05,log(1/1.42870),2),
              c(481,0.02,log(1/1.42870),2),
              c(481,0.01,log(1/1.42870),2),
              c(481,0.1,log(1/1.42870),2) ,
                c(10.58,0.05,log(1.5),1.5),
                c(10.58,0.02,log(1.5),1.5),
                c(10.58,0.01,log(1.5),1.5),
                c(10.58,0.1,log(1.5),1.5))



f1 = function(x){log(exp(2*x^2)-exp(x^2))-(log(37.5^2) -2*481)}
#plot(testd,cont_hill_f(parms,testd))

plot(testd,cont_exp_5_f(Exp5[5,],testd))
plot(testd,cont_exp_5_f(Exp5[6,],testd))
plot(testd,cont_exp_5_f(Exp5[7,],testd))
plot(testd,cont_exp_5_f(Exp5[8,],testd))

mean <- cont_hill_f(parms,doses)
model_list = data.frame(model_list = c(rep("hill",3),rep("exp-3",3),rep("exp-5",3),rep("power",2)),
                        distribution_list =  c(rep(c("normal","normal-ncv","lognormal"),3),"normal",
                                               "normal-ncv"))

y <- rinvgauss(length(mean),mean,4000)

AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                        fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1)
CC <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                        fit_type = "laplace",BMR = 1)

AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                        fit_type = "mcmc",BMD_TYPE = 'hybrid')
CC <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                        fit_type = "laplace",BMD_TYPE = 'hybrid')

AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                        fit_type = "mcmc",BMD_TYPE = 'rel',BMR = 0.1)
CC <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                        fit_type = "laplace",BMD_TYPE = 'rel',BMR = 0.1)

library(splines2)
doses <- seq(0,100,1)
X <- iSpline(doses,knots=seq(30,90,20),intercept=TRUE)
plot(X[,2])
dim(X)

beta <- matrix(rgamma(8,1,2),8,1)
#add +10.58 to get the background
beta1 <- c(0.15486274,0.14054532,0.05806702,0.67421470,0.18371405,0.92821744,1.60669594,1.04508522) 
beta2 <- c(0.677816175,0.322787366,3.356424380,0.102841870,0.009085238,0.087907971,0.344255936,0.100000000)
#add 481 to get the background
beta3 <- c(-45.715480,-53.668821,-14.266726, -8.066441,-1.620326 ,-1.800000,-1.800000 ,-1.800000)
beta4 <- c(-11.9566632,-19.9908719,-39.3450020,-16.3895141,-2.8191037,-4.5801731,-0.834784,-0.7000000)

doses <- seq(0,100,0.5)
X <- iSpline(doses,knots=seq(30,90,20),intercept=TRUE)
bkg = c(10.58,10.58,481,481)
lambda <- c(902.3632,902.3632,78643.17,78643.17)
lambda <- c()
bmd_sd_ig <- c(0,0,0,0)
bmd_sd_no <- c(0,0,0,0)
bmd_sd_ln <- c(0,0,0,0)
i = 1
var_ig <- c(NA,NA,NA,NA)

for(i in 1:4){
  mean = X%*%t(iS[i,,drop=F]) + bkg[i]
  isM <- splinefun(doses,mean)

  var_ig[i] <- isM(0)^3/lambda[i]
  if (isM(0) < isM(100)){
    sdDif = sqrt(var_ig[i])
  }else{
    sdDif = sqrt(var_ig[i])
  }
  
  bb  <-  function(d){abs(isM(d)-isM(0))-sdDif} 
  bmd_sd_ig[i] <- uniroot(bb,c(0,100))$root
}

var_no <- var_ig
for(i in 1:4){
  mean = X%*%t(iS[i,,drop=F]) + bkg[i]
  isM <- splinefun(doses,mean)
 
  if (isM(0) < isM(100)){
    sdDif = sqrt(var_no[i])
  }else{
    sdDif = sqrt(var_no[i])
  }
  
  bb  <-  function(d){abs(isM(d)-isM(0))-sdDif} 
  bmd_sd_no[i] <- uniroot(bb,c(0,100))$root
}

var_ln <- c(0.0777,0.0777,0.107302,0.107302)
for(i in 1:4){
  mean = X%*%t(iS[i,,drop=F]) + bkg[i]
  isM <- splinefun(doses,mean)
  var_m <- expm1(var_ln[i]^2)*exp(2*log(isM(0))+var_ln[i]^2)
  bb  <-  function(d){abs(isM(d)-isM(0))-sqrt(var_m)} 
  bmd_sd_no[i] <- uniroot(bb,c(0,100))$root
}

