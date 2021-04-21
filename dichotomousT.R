library(ToxicR)


o1 <- c(0.1,0.05, -9999)
o2 <- c(1,2)

mData <- matrix(c(0,    0,100,
                  50,   5,100,
                  100, 30,100,
                  150, 65,100,
                  200, 90,100),nrow=5,ncol=3,byrow=T)


mData <- matrix(c(12,   1,5,
                  15,   2,6,
                  18,   3,7,
                  21,   4,8),nrow=4,ncol=3,byrow=T)
Q = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],
                           fit_type = "mcmc")
single_dichotomous_fit()