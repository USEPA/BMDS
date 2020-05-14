library(BMDS)


o1 <- c(0.1,0.05, -9999)
o2 <- c(1,3)

mData <- matrix(c(0,    0,100,
                  50,   5,100,
                  100, 30,100,
                  150, 65,100,
                  200, 90,100),nrow=5,ncol=3,byrow=T)


prior <- matrix(c(1,	0,	2,	-20,	20,
                  1,	log(2), sqrt(0.18)	,	-40,	40,
                  2,	0,	1,	0,	40),nrow=3,ncol=5,byrow=T)
#Ta<-bmd_single_dichotomous("log-probit",mData,o1,o2)
#Tb<-bmd_single_dichotomous("log-probit",mData,o1,o2,PR=prior,isBAYES = F)

bmd_ma_dichotomous(mData,o1,o2)
