library(ToxicR)
library(ggplot2)

library(readr)
PFOA_Liver <- read_table2("~/Documents/PFOA_Liver.txt", 
                         col_names = FALSE)


doses	<- c(0,	0,	0,	0,	0.156,	0.156,	0.156,	0.3125,	0.3125,	0.3125,
           0.625,	0.625,	0.625,	1.25,	1.25,	1.25,	2.5,	2.5,	2.5,	5,5,
           5,	5,	10,	10,	10,	10,	20,	20,	20,	20) 


library(dplyr)
library(ggplot2)
library(ToxicR)

for (ii in 1:5){
  print(ii)

  temp <- PFOA_Liver[ii,]
  v1 <- as.numeric(temp[2:length(temp)])
  library(ToxicR)
  SS<- single_continuous_fit(as.matrix(doses),as.matrix(v1),distribution = "normal-ncv",model_type = "hill",
                              BMD_TYPE = "sd",BMR = 1 , fit_type = "laplace",isFast = F)
   TT <- single_continuous_fit(as.matrix(doses),as.matrix(v1),distribution = "normal-ncv",model_type = "hill",
                            BMD_TYPE = "sd",BMR = 1 , fit_type = "laplace",isFast = T)
  
  RR<- single_continuous_fit(as.matrix(doses),as.matrix(v1),distribution = "normal-ncv",model_type = "hill",
                              BMD_TYPE = "sd",BMR = 1 , fit_type = "mcmc",isFast = F)
  print(BB$bmd)
  plot(BB)
}
