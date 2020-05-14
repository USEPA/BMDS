#set.seed(12345)
library(BMDS) # Uncomment this line if BMDS is not already manually loaded in R console
library(Rcpp)

############################################################################
####################     Data set: Continuous2.dax     #####################
############################################################################

c2ss <- F # Individual dose-response data

c2_doses <- c(0,0,0,0,18,18,18,18,18,20,20,20,20,30,30,30,30,35,35,35,35,40,40,40,40,40)
c2_y <- c(39,38.4,36.3,37.1,40.2,45.3,42.1,38.3,35.9,42.5,45.2,40.1,39.8,50.1,53.4,48.2,52.1,56.1,50.4,53.2,
       55.2,55.1,59.1,56.3,52.9,53.7)
c2_max_dose<-max(c2_doses)

####################     Exponential Model - E3     #####################

# - Normal, 1-SD, CV
c2priorE3_1SD_opt1<-  matrix(c(0, 34.9048, 1, 0, 1e8, # a
                        0, 0.464051, 1, 0, 1e8,     # b
                        0, 4, 1, 0, 1e8,    # log(c)
                        0, 1, 0, 1, 1000, #d 
                        0, 2.19817 , 1, -1000, 1000) # ln-alpha
                      ,5,5,byrow=T)

# - Normal, 1-SD, NCV
c2priorE3_1SD_opt2<-  matrix(c(0, 36.98, 1, 0, 1e8, # a
                         0, 0.584556, 1, 0, 1e8,     # b
                         0, 4, 1, 0, 1e8,    # log(c)
                         0, 1.57314, 1, 1, 1000, #d 
                         0, 0.870443, 1, -1000, 1000, # rho
                         0, -1.36684, 1, -1000, 1000) # ln-alpha
                       ,6,5,byrow=T)  

# - Log-normal, 1-SD 
c2priorE3_1SD_opt3<-  matrix(c(0, 35.4338, 1, 0, 1e8, # a
                        0, 0.437486, 1, 0, 100,     # b
                        0, 7, 1, 0, 100,    # log(c)
                        0, 2, 1, 1, 100, #d 
                        0, -5.33291, 1, -1000, 1000) # ln-alpha
                       ,5,5,byrow=T)

print("Running EXP-3 Opt1")
c2_E3_opt1_res <- bmd_single_continuous('exp-3', 'STDev', BMRF = 1, bkg_prob=-9999,
                                    PR=c2priorE3_1SD_opt1, DATA=cbind(c2_doses,c2_y), sstat=c2ss,
                                    constVar=T, is_log_normal=F, alpha=0.05)

print("EXP-3 Opt1 results:")
print(paste("BMD values:", paste(c2_E3_opt1_res$BMD * c2_max_dose, collapse = ", ")))
print(c2_E3_opt1_res$EST)

#print("Running EXP-3 Opt2")
# c2_E3_opt2_res <- bmd_single_continuous('exp-3', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c2priorE3_1SD_opt2, DATA=cbind(c2_doses,c2_y), sstat=c2ss,
#                                         constVar=F, is_log_normal=F, alpha=0.05)
# 
#print("EXP-3 Opt2 results:")
# print(paste("BMD values:", paste(c2_E3_opt2_res$BMD * c2_max_dose, collapse = ", ")))
# print(c2_E3_opt2_res$EST)

#print("Running EXP-3 Opt3")
# c2_E3_opt3_res <- bmd_single_continuous('exp-3', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c2priorE3_1SD_opt3, DATA=cbind(c2_doses,c2_y), sstat=c2ss,
#                                         constVar=T, is_log_normal=T, alpha=0.05)
# 
#print("EXP-3 Opt3 results:")
# print(paste("BMD values:", paste(c2_E3_opt3_res$BMD * c2_max_dose, collapse = ", ")))
# print(c2_E3_opt3_res$EST)

####################     Exponential Model - E5     #####################

# - Normal, 1-SD, CV
c2priorE5_1SD_opt1<-  matrix(c(0, 34.9048, 1, 0, 1e8, # a
                               0, 0.0021172, 1, 0, 1e8,     # b
                               0, 5.62604, 1, 0, 1e8,    # log(c)
                               0, 1.5301, 0, 1, 1000, #d 
                               0, 2.19817 , 1, -1000, 1000) # ln-alpha
                             ,5,5,byrow=T)

# - Normal, 1-SD, NCV
c2priorE5_1SD_opt2<-  matrix(c(0, 36.98, 1, 0, 1e8, # a
                               0, 0.584556, 1, 0, 1e8,     # b
                               0, 4, 1, 0, 1e8,    # log(c)
                               0, 1.57314, 1, 1, 1000, #d 
                               0, 0.870443, 1, -1000, 1000, # rho
                               0, -1.36684, 1, -1000, 1000) # ln-alpha
                             ,6,5,byrow=T)  

# - Log-normal, 1-SD 
c2priorE5_1SD_opt3<-  matrix(c(0, 35.4338, 1, 0, 1e8, # a
                               0, 0.437486, 1, 0, 100,     # b
                               0, 7, 1, 0, 100,    # log(c)
                               0, 2, 1, 1, 100, #d 
                               0, -5.33291, 1, -1000, 1000) # ln-alpha
                             ,5,5,byrow=T)

print("Running EXP-5 Opt1")
c2_E5_opt1_res <- bmd_single_continuous('exp-5', 'STDev', BMRF = 1, bkg_prob=-9999,
                                        PR=c2priorE5_1SD_opt1, DATA=cbind(c2_doses,c2_y), sstat=c2ss,
                                        constVar=T, is_log_normal=F, alpha=0.05)

print("EXP-5 Opt1 results:")
print(paste("BMD values:", paste(c2_E5_opt1_res$BMD * c2_max_dose, collapse = ", ")))
print(c2_E5_opt1_res$EST)

#print("Running EXP-5 Opt2")
# c2_E5_opt2_res <- bmd_single_continuous('exp-5', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c2priorE5_1SD_opt2, DATA=cbind(c2_doses,c2_y), sstat=c2ss,
#                                         constVar=F, is_log_normal=F, alpha=0.05)
# 
#print("EXP-5 Opt2 results:")
# print(paste("BMD values:", paste(c2_E5_opt2_res$BMD * c2_max_dose, collapse = ", ")))
# print(c2_E5_opt2_res$EST)

#print("Running EXP-5 Opt3")
# c2_E5_opt3_res <- bmd_single_continuous('exp-5', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c2priorE5_1SD_opt3, DATA=cbind(c2_doses,c2_y), sstat=c2ss,
#                                         constVar=T, is_log_normal=T, alpha=0.05)
# 
#print("EXP-5 Opt3 results:")
# print(paste("BMD values:", paste(c2_E5_opt3_res$BMD * c2_max_dose, collapse = ", ")))
# print(c2_E5_opt3_res$EST)

