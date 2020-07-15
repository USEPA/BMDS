library(ToxicR)


o1 <- c(0.1,0.05, -9999)
o2 <- c(1,2)

mData <- matrix(c(0,    0,100,
                  50,   5,100,
                  100, 30,100,
                  150, 65,100,
                  200, 90,100),nrow=5,ncol=3,byrow=T)

system.time({Q = bmd_ma_dichotomous(mData,o1,o2)})
