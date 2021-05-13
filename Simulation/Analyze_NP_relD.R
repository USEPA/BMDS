relD <- c(	64.86	,	64.86	,	64.86	,
           14.27	,	14.27	,	14.27	,
           9.15	,	9.15	,	9.15	,
           27.79	,	27.79	,	27.79	)


BMD.RD <- matrix(relD,nrow=4,ncol=3,byrow = TRUE)

setwd("~/Documents/r_software/RBMDS/Simulation/Non-Parametric/results2")
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

lap1.rd <- array(NA,c(4,3,4))
lap2.rd <- array(NA,c(4,3,4))
mcmc1.rd <- array(NA,c(4,3,4))
mcmc2.rd <- array(NA,c(4,3,4))


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
  true_BMD.rd = BMD.RD[scond,simtype]

  lap1.rd[scond,simtype,even]  = mean(BMD_result_REL_ML1_lapl[,2] < true_BMD.rd,na.rm=TRUE)
  lap2.rd[scond,simtype,even]  = mean(BMD_result_REL_ML2_lapl[,2] < true_BMD.rd,na.rm=TRUE)
  mcmc1.rd[scond,simtype,even] = mean(BMD_result_REL_ML1_mcmc[,2]  < true_BMD.rd,na.rm=TRUE)
  mcmc2.rd[scond,simtype,even] = mean(BMD_result_REL_ML2_mcmc[,2]  < true_BMD.rd,na.rm=TRUE)
  
}

library(xtable)
RD_no_G_4 <- cbind(lap1.rd[,1,3],lap2.rd[,1,3],mcmc1.rd[,1,3],mcmc2.rd[,1,3])*100
RD_no_G_5 <- cbind(lap1.rd[,1,4],lap2.rd[,1,4],mcmc1.rd[,1,4],mcmc2.rd[,1,4])*100

RD_no_E_4 <- cbind(lap1.rd[,1,1],lap2.rd[,1,1],mcmc1.rd[,1,1],mcmc2.rd[,1,1])*100
RD_no_E_5 <- cbind(lap1.rd[,1,2],lap2.rd[,1,2],mcmc1.rd[,1,2],mcmc2.rd[,1,2])*100

RD_lno_G_4 <- cbind(lap1.rd[,2,3],lap2.rd[,2,3],mcmc1.rd[,2,3],mcmc2.rd[,2,3])*100
RD_lno_G_5 <- cbind(lap1.rd[,2,4],lap2.rd[,2,4],mcmc1.rd[,2,4],mcmc2.rd[,2,4])*100

RD_lno_E_4 <- cbind(lap1.rd[,2,1],lap2.rd[,2,1],mcmc1.rd[,2,1],mcmc2.rd[,2,1])*100
RD_lno_E_5 <- cbind(lap1.rd[,2,2],lap2.rd[,2,2],mcmc1.rd[,2,2],mcmc2.rd[,2,2])*100

RD_ig_G_4 <- cbind(lap1.rd[,3,3],lap2.rd[,3,3],mcmc1.rd[,3,3],mcmc2.rd[,3,3])*100
RD_ig_G_5 <- cbind(lap1.rd[,3,4],lap2.rd[,3,4],mcmc1.rd[,3,4],mcmc2.rd[,3,4])*100

RD_ig_E_4 <- cbind(lap1.rd[,3,1],lap2.rd[,3,1],mcmc1.rd[,3,1],mcmc2.rd[,3,1])*100
RD_ig_E_5 <- cbind(lap1.rd[,3,2],lap2.rd[,3,2],mcmc1.rd[,3,2],mcmc2.rd[,3,2])*100

xtable(cbind(RD_no_E_4,RD_no_G_4), digits = 1)
xtable(cbind(RD_no_E_5,RD_no_G_5), digits = 1)

xtable(cbind(RD_lno_E_4,RD_lno_G_4), digits = 1)
xtable(cbind(RD_lno_E_5,RD_lno_G_5), digits = 1)

xtable(cbind(RD_ig_E_4,RD_ig_G_4), digits = 1)
xtable(cbind(RD_ig_E_5,RD_ig_G_5), digits = 1)

