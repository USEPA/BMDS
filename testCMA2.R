#library(ggplot2)
#library(reshape)
#library(cowplot)
library(Rcpp)
library(BMDS)
############################################################################
####################     Data set: ??????     #####################
############################################################################

c2ss <- F # Individual dose-response data
c2_doses <- c(0 , 0 , 0 , 0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,
              37.5 ,  37.5 ,37.5 ,37.5 ,37.5 ,37.5 ,37.5 ,  37.5 ,
              37.5 ,  37.5 ,75 ,75 ,75 ,75 ,75 ,75 ,75 ,75 ,  75 ,
              75 ,150 ,150,150,150,150,150,150,150,150,150 ,
              300 ,300,300,300,300,300,300,300, 300 , 300 ,
              600 , 600 , 600 , 600 , 600 , 600, 600, 600 ,
              600 , 600 );

c2_y<-c(3.992,3.766,4.071,4.096,4.047,3.745,3.822,4.580,3.836,3.997,4.292,4.289,4.473,4.518,4.176,3.895,4.009,3.809,4.181,4.073,
        4.161,4.009,3.968,3.888,3.905,4.495,4.194,3.949,3.994,4.005,4.559,4.325,4.286,4.521,
        4.711,4.013,4.568,4.216,4.292,4.118,5.660,5.995,6.510,6.756,6.517,7.566,6.176,6.826,
        7.429,6.769,10.403,11.745,9.977,11.569,11.334,12.586,10.007,12.011,10.388,12.245)

fit <- bmd_ma_con2('STDev',cbind(c2_doses,c2_y),sstat=F)
mydata<-data.frame(c2_doses,c2_y)
class(fit) <- "doseResponse"
#plot_ma_continuous(fit)
class(fit) <- "bmdCDF"
#plot_ma_continuous(fit) #.bmdCDF(fit)
class(fit) <- "summaryTable"
print_ma_continuous(fit)
print("StDev results:")
#sprintf("BMD:%3.2f (%3.2f,%3.2f)",get_ma_bmd(fit,0.5),get_ma_bmd(fit,0.05),get_ma_bmd(fit,0.95))
print(paste0("BMD=",get_ma_bmd(fit,0.5)," (",get_ma_bmd(fit,0.05),", ",get_ma_bmd(fit,0.95), ")"))

fit <- bmd_ma_con2('RelDev',cbind(c2_doses,c2_y),sstat=F,BMRF=0.05)
#fit <- bmd_ma_continuous('RelDev',cbind(c2_doses,c2_y),sstat=F,BMRF=0.05)
mydata<-data.frame(c2_doses,c2_y)
class(fit) <- "doseResponse"
#plot_ma_continuous(fit)
class(fit) <- "bmdCDF"
#plot_ma_continuous.bmdCDF(fit)
class(fit) <- "summaryTable"
print_ma_continuous(fit)
print("RelDev results:")
#sprintf("BMD:%3.2f (%3.2f,%3.2f)",get_ma_bmd(fit,0.5),get_ma_bmd(fit,0.05),get_ma_bmd(fit,0.95))
print(paste0("BMD=",get_ma_bmd(fit,0.5)," (",get_ma_bmd(fit,0.05),", ",get_ma_bmd(fit,0.95), ")"))

# BMD <- fit$ma$MA_BMD_CDF
# BMD[,1] = BMD[,1]*fit$dose_scale
# bmd_ma = splinefun(BMD[,2],BMD[,1])
# X <- bmd_ma(runif(100000))
# hist(X,main="BMD Rel-Dev 5%",xlab="BMD")
