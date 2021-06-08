
#06/07/21 SL update 
library(ToxicR)
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

# Adjust size of the interval bar;
# Changed error bar with reasonable size 
# However, there is still some error bar gap 
# I think it's inevitable have some space between 

plot(c)
.plot.BMDcont_fit_maximized(c,qprob=0.05)



#MCMC test 
c2 = single_continuous_fit(M2[,1,drop=F],M2[,2:4],sstat=F,BMD_TYPE="sd",BMR=1, 
                          distribution = "normal",fit_type="mcmc",model_type = "power")


# 06/07/21 SL -- tested and updated for MCMC - single continous fit 

plot(c2)
.plot.BMDcont_fit_MCMC(c2,qprob=0.05)


# Data generation for testing purpose


model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))



# 05/28 SL Try to focus on the issue here
# 1. CI Band seems odd
# 2. BMD Density plot looks odd 
# 3. Color Theme update for giving user more option

# I think it is useful that which model is dominant in terms of the line- color
# Range max-min part should be updated... 
# BMDU is out of bound


# How about the log scale?
c3 <- ma_continuous_fit(M2[,1,drop=F],M2[,2:4],fit_type = "mcmc",BMR = 2,model_list = model_list )

plot(c3)
.plot.BMDcontinuous_MA(c3)






# This is not a sufficent statistics one. BMDL ? is way too high? I guess
# If the test doses are high and BMD


doses = c( 0, 0, 0, 0, 18, 18, 18, 18, 18, 20, 20, 20, 20, 30, 30, 30, 30, 35, 35, 35, 35, 40, 40, 40, 40, 40);
Y = c(39,38.4,36.3,37.1, 40.2,45.3,42.1,38.3,35.9, 42.5,45.2,40.1,39.8, 50.1,53.4,48.2,52.1, 56.1,50.4,53.2,55.2, 35.1,39.1,36.3,32.9,33.7);
Y[doses==20] = Y[doses==20] + 5.6
B  <- single_continuous_fit(as.matrix(doses),as.matrix(Y),model_type = "polynomial", distribution="normal",fit_type = "laplace",degree=3,BMR = 1,sstat = F,samples = 150000)

# SL - This is not a sufficient Statistics so we don't need to add 
plot(B)
.plot.BMDcont_fit_maximized(B,qprob=0.05)
