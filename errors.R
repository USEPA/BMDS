#06/07/21 SL update 
library(ToxicR)
M2           <- matrix(0,nrow=5,ncol=4)

# double D[] = {0,50, 100, 150, 200};
# double Y[] = {10, 0 , -10, -20 ,-30};
# double N[] = {100, 100, 100, 100, 100};
# double SD[] = {3, 4, 5, 6, 7};
colnames(M2) <- c("Dose","Resp","N","StDev")
#M2[, 1]      <- c(0,	0.156,	0.312,	0.625,	1.25,	2.5)
#M2[, 2]      <- c(33.52,	37.66,	40.08,	44.25,	50.84,	67.75)
#M2[, 3]      <- c(10,	10,	10,	10,	10,	10)
#M2[, 4]      <- c(2.371708245,	2.814427118,	1.77087549,	2.593067681,	2.118726032,	2.846049894)

M2[,1] <- c(0,25,50,100,200)
M2[,3] <- c(20,20,19,20,20)
M2[,2] <- c(6,5.2,2.4,1.1,0.75)
M2[,4]<-  c(1.2,1.1,0.81,0.74,0.66)
#double D[] = {0,50, 100, 150, 200};
#double Y[] = {10, 20 , 30, 40 ,50};
#double N[] = {100, 100, 100, 100, 100};
#double SD[] = {3, 4, 5, 6, 7};

#fast ewald estimate
c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                          distribution = "normal",fit_type="laplace",model_type = "hill",degree = 4,ewald=F)
plot(c)


#MCMC estimate 
c2 = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1,
                          distribution = "normal",fit_type="mcmc",model_type = "exp-5")

plot(c2)

cleveland_plot(c2)
BMDcontinous_MA_MCMC(c2)






