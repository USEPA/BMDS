library(ToxicR)
library(readr)

doses<- c(0,0,0,0,0.156,0.156,0.156,0.3125,0.3125,0.3125,
          0.625,0.625,0.625,1.25,1.25,1.25,2.5,2.5,2.5,5,5,
          5,5,10,10,10,10,20,20,20,20) 

PFOA_Liver <- read_table2("PFOA_Liver.txt", 
                          col_names = FALSE)

library(dplyr)

temp <- PFOA_Liver %>% filter(X1 == "CYP3A1_32809")
v1 <- as.numeric(temp[2:length(temp)])
fit1  <- single_continuous_fit(as.matrix((doses)),as.matrix(v1),model_type = "exp-5", distribution="normal",degree=3,fit_type = "laplace",BMR = 3,isFast=T)
fit2  <- single_continuous_fit(as.matrix((doses)),as.matrix(v1),model_type = "hill", distribution="normal",
                               degree=3,fit_type = "mcmc",BMR = 0.05,BMD_TYPE = "rel") #5% increase 
fit3  <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "laplace",BMD_TYPE = "sd", BMR =1.5) #1.5 standard deviations



