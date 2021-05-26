#set.seed(12345)
#library(BMDS) # Uncomment this line if BMDS is not already manually loaded in R console
#library(Rcpp)

############################################################################
####################     Data set: Continuous3.dax     #####################
############################################################################

c3ss <- T # Summarized dose-response data

c3Dat <- matrix(0,nrow=5,ncol=4)
colnames(c3Dat) <- c("Dose","Resp","N","StDev")
c3Dat[, 1] <- c(0, 35, 105, 316, 625)
c3Dat[, 2] <- c(1.61, 1.66, 1.75, 1.81, 1.89)
c3Dat[, 3] <- 2
c3Dat[, 4] <- c(0.12, 0.13, 0.11, 0.15, 0.13)
c3Dat[, 4] <- 5
c3Dat[, 3] <- c(0.12, 0.13, 0.11, 0.15, 0.13)
c3_max_dose<-max(c3Dat[, 1])


c3Dat


print("Running EXP-3 Opt1")
BB <- single_continuous_fit(as.matrix(c3Dat[, 1]),c3Dat[, 2:4],model_type = "exp-5",
                            distribution="normal-ncv",fit_type = "mle",sstat = T)


print("EXP-3 Opt1 results:")
paste("BMD values:", paste(c3_E3_opt1_res$BMD * c3_max_dose, collapse = ", "))
c3_E3_opt1_res$EST

#print("Running EXP-3 Opt2")
# c3_E3_opt2_res <- bmd_single_continuous('exp-3', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c3priorE3_1SD_opt2, DATA=c3Dat, sstat=c3ss,
#                                         constVar=F, is_log_normal=F, alpha=0.05)
# 
#print("EXP-3 Opt2 results:")
# paste("BMD values:", paste(c3_E3_opt2_res$BMD * c3_max_dose, collapse = ", "))
# c3_E3_opt2_res$EST

#print("Running EXP-3 Opt3")
# c3_E3_opt3_res <- bmd_single_continuous('exp-3', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c3priorE3_1SD_opt3, DATA=c3Dat, sstat=c3ss,
#                                         constVar=T, is_log_normal=T, alpha=0.05)
# 
#print("EXP-3 Opt3 results:")
# paste("BMD values:", paste(c3_E3_opt3_res$BMD * c3_max_dose, collapse = ", "))
# c3_E3_opt3_res$EST

####################     Exponential Model - E5     #####################

# - Normal, 1-SD, CV
c3priorE5_1SD_opt1<-    matrix(c(0,34.9048,1, 0,1e8, # a
                                 0,0.464049,1, 0,1e8,     # b
                                 0,4,1, 0,1e8,    # log(c)
                                 0,1,0,1,1000, #d 
                                 0,2.19817,1,-1000,1000) # ln-alpha
                               ,5,5,byrow=T)

# - Normal, 1-SD, NCV
c3priorE5_1SD_opt2<-  matrix(c(0,36.98,1, 0,1e8, # a
                               0,0.584556,1, 0,1e8,     # b
                               0,4,1, 0,1e8,    # log(c)
                               0,1.57314,1,1,1000, #d 
                               0,0.870443,1,-1000,1000, # rho
                               0,-1.36684,1,-1000,1000) # ln-alpha
                             ,6,5,byrow=T)  

# - Log-normal, 1-SD 
c3priorE5_1SD_opt3<-    matrix(c(0,35.4338,1, 0,1e8, # a
                                 0,0.437486,1, 0,100,     # b
                                 0,7,1, 0,100,    # log(c)
                                 0,2,1,1,100, #d 
                                 0,-5.33291,1,-1000,1000) # ln-alpha
                               ,5,5,byrow=T)


print("Running EXP-5 Opt1")
c3_E5_opt1_res <- bmd_single_continuous('exp-5', 'STDev', BMRF = 1, bkg_prob=-9999,
                                        PR=c3priorE5_1SD_opt1, DATA=c3Dat, sstat=c3ss,
                                        constVar=T, is_log_normal=F, alpha=0.05)

print("EXP-5 Opt1 results:")
paste("BMD values:", paste(c3_E5_opt1_res$BMD * c3_max_dose, collapse = ", "))
c3_E5_opt1_res$EST

#print("Running EXP-5 Opt2")
# c3_E5_opt2_res <- bmd_single_continuous('exp-5', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c3priorE5_1SD_opt2, DATA=c3Dat, sstat=c3ss,
#                                         constVar=F, is_log_normal=F, alpha=0.05)
# 
#print("EXP-5 Opt2 results:")
# paste("BMD values:", paste(c3_E5_opt2_res$BMD * c3_max_dose, collapse = ", "))
# c3_E5_opt2_res$EST

#print("Running EXP-5 Opt3")
# c3_E5_opt3_res <- bmd_single_continuous('exp-5', 'STDev', BMRF = 1, bkg_prob=-9999,
#                                         PR=c3priorE5_1SD_opt3, DATA=c3Dat, sstat=c3ss,
#                                         constVar=T, is_log_normal=T, alpha=0.05)
# 
#print("EXP-5 Opt3 results:")
# paste("BMD values:", paste(c3_E5_opt3_res$BMD * c3_max_dose, collapse = ", "))
# c3_E5_opt3_res$EST

