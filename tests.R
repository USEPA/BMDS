library(ToxicR)

library(readr)
library(dplyr)
PFOA_Liver <- read_table2("~/Documents/PFOA_Liver.txt", 
                          col_names = FALSE)

doses <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.948669,0.948669,0.948669,2.999955,2.999955,2.999955,9.48669,9.48669,9.48669,29.999548,29.999548,29.999548,94.8669,94.8669,94.8669,299.99548,299.99548,299.99548,948.669,948.669,948.669,2999.9548,2999.9548,2999.9548,9486.69,9486.69,9486.69)
v2 <- c(6.957,7.838,7.157,7.589,7.808,7.305,7.143,7.778,7.449,7.761,7.485,8.03,6.863,7.062,7.611,7.135,7.484,7.161,7.272,7.603,7.428,7.291,7.048,7.562,7.324,7.17,7.072,7.252,7.249,7.397,6.983,7.175,7.13,8.341,8.289,8.201,8.838,8.865,8.971)


R  <- single_continuous_fit(as.matrix(doses),as.matrix(v2),model_type = "polynomial", 
                            distribution="normal",degree=4,fit_type = "laplace",BMR = 3,isFast = T)

R  <- single_continuous_fit(as.matrix(doses),as.matrix(v2),model_type = "power", 
                            distribution="normal",degree=4,fit_type = "laplace",BMR = 3,isFast = T)

R  <- single_continuous_fit(as.matrix(doses),as.matrix(v2),model_type = "exp-5", 
                            distribution="normal-ncv",degree=4,fit_type = "laplace",BMR = 3,isFast = T)


R  <- single_continuous_fit(as.matrix(doses),as.matrix(v2),model_type = "hill", 
                            distribution="normal",degree=4,fit_type = "laplace",BMR = 3,isFast = T)

R  <- single_continuous_fit(as.matrix(doses),as.matrix(v2),model_type = "exp-5", 
                            distribution="lognormal",degree=4,fit_type = "laplace",BMR = 3,isFast = T)


BB <- ma_continuous_fit(as.matrix(doses),as.matrix(v2),fit_type = "laplace",BMR = 1)



S  <- single_continuous_fit(as.matrix(doses),as.matrix(v2),model_type = "hill", 
                            distribution="normal",degree=4,fit_type = "mle",BMR = 3)
