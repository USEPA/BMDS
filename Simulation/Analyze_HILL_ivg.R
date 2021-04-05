bmds.SD<- c(54.9,	54.91,	54.9,
            21.59,	21.6,	21.59,
            7.24,	7.25,	7.24,
            40.92,	40.93,	40.92,
            49.97,	49.97,	49.97,
            16.87,	16.87,	16.87,
            8.31,	8.31,	8.31,
            37.23,	37.23,	37.23)

bmds.H <- c(41.49,	57.19,	42.61,
            10.61,	23.94,	11.35,
            3.13,	8.18,	3.39,
            32.48,	42.32,	33.2,
            40.08,	54.81,	38.9,
            13.04,	18.79,	12.6,
            5.65,	9.78,	5.37,
            30.7,	40.37,	29.91)

BMD.H  <- matrix(bmds.H,nrow=8,ncol=3,byrow=T)
BMD.SD <- matrix(bmds.SD,nrow=8,ncol=3,byrow=T)
setwd("~/Documents/r_software/RBMDS/Simulation/Hill/results")
files <- dir()

ivgSims = grepl("invGaussian",files)
norSims = grepl("normal",files)
lnorSims = grepl("lognormal",files)

cond1 = grepl("sim_1",files)
cond2 = grepl("sim_2",files)
cond3 = grepl("sim_3",files)
cond4 = grepl("sim_4",files)
cond5 = grepl("sim_5",files)
cond6 = grepl("sim_6",files)
cond7 = grepl("sim_7",files)
cond8 = !(cond1 | cond2 | cond3 | cond4 | cond5 | cond6 | cond7)

is_geom = grepl("_g_",files)
is_4    = grepl("_4_",files)

lap1.sd <- array(NA,c(8,3,4))
lap2.sd <- array(NA,c(8,3,4))
mcmc1.sd <- array(NA,c(8,3,4))
mcmc2.sd <- array(NA,c(8,3,4))

lap1.h <- array(NA,c(8,3,4))
lap2.h <- array(NA,c(8,3,4))
mcmc1.h <- array(NA,c(8,3,4))
mcmc2.h <- array(NA,c(8,3,4))

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
  scond = which(c(cond1[ii],cond2[ii],cond3[ii],cond4[ii],cond5[ii],cond6[ii],cond7[ii],cond8[ii]))
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
SD_no_G_4 <- cbind(lap1.sd[,1,3],lap2.sd[,1,3],mcmc1.sd[,1,3],mcmc2.sd[,1,3])*100
SD_no_G_5 <- cbind(lap1.sd[,1,4],lap2.sd[,1,4],mcmc1.sd[,1,4],mcmc2.sd[,1,4])*100

SD_no_E_4 <- cbind(lap1.sd[,1,1],lap2.sd[,1,1],mcmc1.sd[,1,1],mcmc2.sd[,1,1])*100
SD_no_E_5 <- cbind(lap1.sd[,1,2],lap2.sd[,1,2],mcmc1.sd[,1,2],mcmc2.sd[,1,2])*100

SD_lno_G_4 <- cbind(lap1.sd[,2,3],lap2.sd[,2,3],mcmc1.sd[,2,3],mcmc2.sd[,2,3])*100
SD_lno_G_5 <- cbind(lap1.sd[,2,4],lap2.sd[,2,4],mcmc1.sd[,2,4],mcmc2.sd[,2,4])*100

SD_lno_E_4 <- cbind(lap1.sd[,2,1],lap2.sd[,2,1],mcmc1.sd[,2,1],mcmc2.sd[,2,1])*100
SD_lno_E_5 <- cbind(lap1.sd[,2,2],lap2.sd[,2,2],mcmc1.sd[,2,2],mcmc2.sd[,2,2])*100

SD_ig_G_4 <- cbind(lap1.sd[,3,3],lap2.sd[,3,3],mcmc1.sd[,3,3],mcmc2.sd[,3,3])*100
SD_ig_G_5 <- cbind(lap1.sd[,3,4],lap2.sd[,3,4],mcmc1.sd[,3,4],mcmc2.sd[,3,4])*100

SD_ig_E_4 <- cbind(lap1.sd[,3,1],lap2.sd[,3,1],mcmc1.sd[,3,1],mcmc2.sd[,3,1])*100
SD_ig_E_5 <- cbind(lap1.sd[,3,2],lap2.sd[,3,2],mcmc1.sd[,3,2],mcmc2.sd[,3,2])*100

xtable(cbind(SD_no_E_4,SD_no_G_4), digits = 1)
xtable(cbind(SD_no_E_5,SD_no_G_5), digits = 1)

xtable(cbind(SD_lno_E_4,SD_lno_G_4), digits = 1)
xtable(cbind(SD_lno_E_5,SD_lno_G_5), digits = 1)

xtable(cbind(SD_ig_E_4,SD_ig_G_4), digits = 1)
xtable(cbind(SD_ig_E_5,SD_ig_G_5), digits = 1)


HB_ig_G_4 <- cbind(lap1.h[,3,3],lap2.h[,3,3],mcmc1.h[,3,3],mcmc2.h[,3,3])*100
HB_ig_G_5 <- cbind(lap1.h[,3,4],lap2.h[,3,4],mcmc1.h[,3,4],mcmc2.h[,3,4])*100

HB_ig_E_4 <- cbind(lap1.h[,3,1],lap2.h[,3,1],mcmc1.h[,3,1],mcmc2.h[,3,1])*100
HB_ig_E_5 <- cbind(lap1.h[,3,2],lap2.h[,3,2],mcmc1.h[,3,2],mcmc2.h[,3,2])*100

HB_lno_G_4 <- cbind(lap1.h[,2,3],lap2.h[,2,3],mcmc1.h[,2,3],mcmc2.h[,2,3])*100
HB_lno_G_5 <- cbind(lap1.h[,2,4],lap2.h[,2,4],mcmc1.h[,2,4],mcmc2.h[,2,4])*100

HB_lno_E_4 <- cbind(lap1.h[,2,1],lap2.h[,2,1],mcmc1.h[,2,1],mcmc2.h[,2,1])*100
HB_lno_E_5 <- cbind(lap1.h[,2,2],lap2.h[,2,2],mcmc1.h[,2,2],mcmc2.h[,2,2])*100

HB_no_G_4 <- cbind(lap1.h[,1,3],lap2.h[,1,3],mcmc1.h[,1,3],mcmc2.h[,1,3])*100
HB_no_G_5 <- cbind(lap1.h[,1,4],lap2.h[,1,4],mcmc1.h[,1,4],mcmc2.h[,1,4])*100

HB_no_E_4 <- cbind(lap1.h[,1,1],lap2.h[,1,1],mcmc1.h[,1,1],mcmc2.h[,1,1])*100
HB_no_E_5 <- cbind(lap1.h[,1,2],lap2.h[,1,2],mcmc1.h[,1,2],mcmc2.h[,1,2])*100


xtable(cbind(HB_no_E_4,HB_no_G_4), digits = 1)
xtable(cbind(HB_no_E_5,HB_no_G_5), digits = 1)

xtable(cbind(HB_lno_E_4,HB_lno_G_4), digits = 1)
xtable(cbind(HB_lno_E_5,HB_lno_G_5), digits = 1)

xtable(cbind(HB_ig_E_4,HB_ig_G_4), digits = 1)
xtable(cbind(HB_ig_E_5,HB_ig_G_5), digits = 1)

