data90 <- c(0,	    10,	352.2,	19.9,
            125,	  10,	350.6,	11.4,
            250	,   10,	338.8,	20.3,
            500,	  10,	343.5,	15.2,
            1000,	  10,	330.1,	25,
            1500,	  10,	312.5,	21.6)

AA <- matrix(data90,ncol=4,nrow=6,byrow=T)

Y = matrix(NA,nrow=6,ncol=3)
Y[,1] = AA[,3]
Y[,2] = AA[,2]
Y[,3] = AA[,4]
doses<- AA[,1,drop=F]

#double D[] = {0,50, 100, 150, 200};
#double Y[] = {10, 20 , 30, 40 ,50};
#double N[] = {100, 100, 100, 100, 100};
#double SD[] = {3, 4, 5, 6, 7};
library(ToxicR)

c = single_continuous_fit(doses,Y,BMD_TYPE="sd",BMR=1, 
                          distribution = "normal",fit_type= "mle" ,model_type = "polynomial",degree=5)

