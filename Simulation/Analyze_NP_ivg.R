bmds.H <- c(70.97,	72.03,	66.76,
          16.06,	16.41,	14.8,
           6.93,	6.49,	7.54,
          22.15,	20.99,	23.74)

bmds.SD <- c(88.55,	88.55,	88.55,
             26.5,	26.5,	26.5,
             17.12,	17.14,	17.12,
             48.49,	48.54,	48.49)

BMD.H  <- matrix(bmds.H,nrow=4,ncol=3,byrow=T)
BMD.SD <- matrix(bmds.SD,nrow=4,ncol=3,byrow=T)
setwd("~/Documents/r_software/RBMDS/Simulation/Non-Parametric/results")
files <- dir()

ivgSims = grepl("invGaussian",files)
norSims = grepl("normal",files)
lnorSims = grepl("lognormal",files)

cond1 = grepl("sim_1",files)
cond2 = grepl("sim_2",files)
cond3 = grepl("sim_3",files)
cond4 = !(cond1 | cond2 | cond3)
is_geom = grepl("_g_",files)
is_4    = grepl("_4_",files)

lap1.sd <- array(NA,c(4,3,4))
lap2.sd <- array(NA,c(4,3,4))
mcmc1.sd <- array(NA,c(4,3,4))
mcmc2.sd <- array(NA,c(4,3,4))

lap1.h <- array(NA,c(4,3,4))
lap2.h <- array(NA,c(4,3,4))
mcmc1.h <- array(NA,c(4,3,4))
mcmc2.h <- array(NA,c(4,3,4))

for (ii in 1:length(files)){
  load(files[ii]) 
  simtype = 1
  if (lnorSims[ii] == TRUE){
    simtype = 2
  }
  if (ivgSims[ii] == TRUE){
    simtype = 3
  }
  
  even = 1
  if (is_geom[ii] == TRUE){
    even = 3
  }
  if (!is_4[ii] == TRUE){ # for 5 dose groups
    even = even + 1
  }
  #simulation condition
  scond = which(c(cond1[ii],cond2[ii],cond3[ii],cond4[ii]))
  true_BMD.h  = BMD.H[scond,simtype]
  true_BMD.sd = BMD.SD[scond,simtype]
  
  lap1.sd[scond,simtype,even]  = mean(BMD_result_SD_ML1_lapl[,2] < true_BMD.sd,na.rm=TRUE)
  lap2.sd[scond,simtype,even]  = mean(BMD_result_SD_ML2_lapl[,2] < true_BMD.sd,na.rm=TRUE)
  mcmc1.sd[scond,simtype,even] = mean(BMD_result_SD_ML1_mcmc[,2]  < true_BMD.sd,na.rm=TRUE)
  mcmc2.sd[scond,simtype,even] = mean(BMD_result_SD_ML2_mcmc[,2]  < true_BMD.sd,na.rm=TRUE)
  
  lap1.h[scond,simtype,even]   = mean(BMD_result_HB_ML1_lapl[,2] < true_BMD.h,na.rm=TRUE)
  lap2.h[scond,simtype,even]   = mean(BMD_result_HB_ML2_lapl[,2] < true_BMD.h,na.rm=TRUE)
  mcmc1.h[scond,simtype,even]  = mean(BMD_result_HB_ML1_mcmc[,2]  < true_BMD.h,na.rm=TRUE)
  mcmc2.h[scond,simtype,even]  = mean(BMD_result_HB_ML2_mcmc[,2]  < true_BMD.h,na.rm=TRUE)
  
}


library(xtable)

SD_ig_G_4 <- cbind(lap1.sd[,3,3],lap2.sd[,3,3],mcmc1.sd[,3,3],mcmc2.sd[,3,3])*100
SD_ig_G_5 <- cbind(lap1.sd[,3,4],lap2.sd[,3,4],mcmc1.sd[,3,4],mcmc2.sd[,3,4])*100

SD_ig_E_4 <- cbind(lap1.sd[,3,1],lap2.sd[,3,1],mcmc1.sd[,3,1],mcmc2.sd[,3,1])*100
SD_ig_E_5 <- cbind(lap1.sd[,3,2],lap2.sd[,3,2],mcmc1.sd[,3,2],mcmc2.sd[,3,2])*100

xtable(cbind(SD_ig_E_5,SD_ig_G_5), digits = 1)

HB_ig_G_4 <- cbind(lap1.h[,3,3],lap2.h[,3,3],mcmc1.h[,3,3],mcmc2.h[,3,3])*100
HB_ig_G_5 <- cbind(lap1.h[,3,4],lap2.h[,3,4],mcmc1.h[,3,4],mcmc2.h[,3,4])*100

HB_ig_E_4 <- cbind(lap1.h[,3,1],lap2.h[,3,1],mcmc1.h[,3,1],mcmc2.h[,3,1])*100
HB_ig_E_5 <- cbind(lap1.h[,3,2],lap2.h[,3,2],mcmc1.h[,3,2],mcmc2.h[,3,2])*100


xtable(cbind(HB_ig_E_5,HB_ig_G_5), digits = 1)

