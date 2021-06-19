# Updated contents from previous meeting


# 1. Fixed pointy density plot

library(dplyr)
library(ToxicR)
library(readr)
PFOA_Liver <- read_table2("PFOA_Liver.txt", 
                          col_names = FALSE)

temp <- PFOA_Liver %>% filter(X1 == "ABCG2_32656")
v1 <- as.numeric(temp[2:length(temp)])


doses	<- c(0,	0,	0,	0,	0.156,	0.156,	0.156,	0.3125,	0.3125,	0.3125,
           0.625,	0.625,	0.625,	1.25,	1.25,	1.25,	2.5,	2.5,	2.5,	5,5,
           5,	5,	10,	10,	10,	10,	20,	20,	20,	20) 


model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))


BB <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "mcmc",BMR = 2,model_list = model_list )

plot(BB)

# 06/19/21 Test

# Need to distinguish the names for each models
.cleveland_plot.BMDcontinous_MA(BB)


# Question not sure why BMD is superimposed? 
.plot.BMDcontinuous_MA(BB, qprob=0.05)


# 2. Error bar fix for sufficient statistics - continuous

M2           <- matrix(0,nrow=5,ncol=4)
colnames(M2) <- c("Dose","Resp","N","StDev")
M2[, 1]      <- c(0,50, 100, 150, 200)
M2[, 2]      <- c(10, 20 , 30, 40 ,50)
M2[, 3]      <- c(100, 100, 100, 100, 100)
M2[, 4]      <- c(3, 4, 5, 6, 7)
#double D[] = {0,50, 100, 150, 200};
#double Y[] = {10, 20 , 30, 40 ,50};
#double N[] = {100, 100, 100, 100, 100};
#double SD[] = {3, 4, 5, 6, 7};
c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],sstat=F,BMD_TYPE="sd",BMR=1, 
                          distribution = "normal",fit_type="mle",model_type = "power")


plot(c)
.plot.BMDcont_fit_maximized(c,qprob=0.05)



#MCMC test 
c2 = single_continuous_fit(M2[,1,drop=F],M2[,2:4],sstat=F,BMD_TYPE="sd",BMR=1, 
                           distribution = "normal",fit_type="mcmc",model_type = "power")


# 06/07/21 SL -- tested and updated for MCMC - single continuous fit 

plot(c2)
.plot.BMDcont_fit_MCMC(c2,qprob=0.05)


# MA continuous fit -- before the function is not wroking but now it is fine 


model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))

c3 <- ma_continuous_fit(M2[,1,drop=F],M2[,2:4],fit_type = "mcmc",BMR = 2,model_list = model_list )
plot(c3)

# Fix n=512
.plot.BMDcontinuous_MA(c3)



# 3. Dichotomous density color matching with continous density plot 





# Sample Data - Dichotomous Example
mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)


# Single model fitting- MCMC
# too much white space- this needs to be adjusted automatically - Not the for the first Priority, but needs to be updated once 
# Continous part is updated

A_single_mcmc<-single_dichotomous_fit(mData[,1],mData[,2],mData[,3], model_type="hill",fit_type="mcmc")

#Need to color Match with Continous plot 
plot(A_single_mcmc)
.plot.BMDdich_fit_MCMC(A_single_mcmc)+scale_y
