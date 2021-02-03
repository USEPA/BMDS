library(actuar)
library(ToxicR)



x <- rinvgauss(1000, 10, shape = 200)

doses <- rep(c(0,25,50,75,100),each=5)
testd <- seq(0,100,1)

parms <- c(4,5.02,70,3.3)
parms <- c(4,5.02,40,1.3)
parms <- c(10.88,-3.02,25,3)


plot(testd,cont_hill_f(parms,testd))

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