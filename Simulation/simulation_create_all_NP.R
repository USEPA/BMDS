#Create MA simulation
#Simulation Creation for the CMA paper
#iSpline copied from Simulation Conditions.Rmd (2/5/11)
library(ToxicR)
library(actuar)

set.seed(90101)
setwd('./Non-Parametric')
#log-transformed for fitting as in ToxicR use
beta1 <- c(0.15486274,0.14054532,0.05806702,0.67421470,0.18371405,0.92821744,1.60669594,1.04508522) 
beta2 <- c(0.677816175,0.322787366,3.356424380,0.102841870,0.009085238,0.087907971,0.344255936,0.100000000)
#add 481 to get the background
beta3 <- c(-45.715480,-53.668821,-14.266726, -8.066441,-1.620326 ,-1.800000,-1.800000 ,-1.800000)
beta4 <- c(-11.9566632,-19.9908719,-39.3450020,-16.3895141,-2.8191037,-4.5801731,-0.834784,-0.7000000)
iS <- rbind(beta1,beta2,beta3,beta4)

dose_g_5 <- rep(c(0,6.25,12.5,25,50,100),each=10) 
dose_g_5 <- rep(c(0,6.25,12.5,25,50,100),each=10) 
dose_e_5 <- rep(c(0,20,40,60,80,100),each=10) 
dose_g_4 <- rep(c(0,12.5,25,50,100),each=10) 
dose_e_4 <- rep(c(0,25,50,75,100),each=10) 

library(splines2)
X_g_5 <- iSpline(dose_g_5,knots=seq(30,90,20))
X_e_5 <- iSpline(dose_e_5,knots=seq(30,90,20))
X_g_4 <- iSpline(dose_g_4,knots=seq(30,90,20))
X_e_4 <- iSpline(dose_e_4,knots=seq(30,90,20))
bkg = c(10.58,10.58,481,481)
sd    <- c(rep(2.28,2),rep(77.5,2))
#normal data
for (i in 1:4){
  means_g_5 <- X_g_5%*%t(iS[i,,drop=F]) + bkg[i]
  means_e_5 <- X_e_5%*%t(iS[i,,drop=F]) + bkg[i]
  means_g_4 <- X_g_4%*%t(iS[i,,drop=F]) + bkg[i]
  means_e_4 <- X_e_4%*%t(iS[i,,drop=F]) + bkg[i]
  
  doses = dose_g_5
  sim_data <- matrix(rnorm(1000*length(dose_g_5),means_g_5,sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("NP_normal_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rnorm(1000*length(dose_e_5),means_e_5,sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("NP_normal_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rnorm(1000*length(dose_g_4),means_g_4,sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("NP_normal_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rnorm(1000*length(dose_e_4),means_e_4,sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("NP_normal_e_4_sim_%d.Rdata",i))
}
#log-normal data
sd    <- c(rep(0.209,2),rep(0.158,2))
for (i in 1:4){
  means_g_5 <- X_g_5%*%t(iS[i,,drop=F]) + bkg[i]
  means_e_5 <- X_e_5%*%t(iS[i,,drop=F]) + bkg[i]
  means_g_4 <- X_g_4%*%t(iS[i,,drop=F]) + bkg[i]
  means_e_4 <- X_e_4%*%t(iS[i,,drop=F]) + bkg[i]
  
  
  doses = dose_g_5
  sim_data <- matrix(rlnorm(1000*length(dose_g_5),log(means_g_5),sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("NP_lognormal_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rlnorm(1000*length(dose_e_5),log(means_e_5),sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("NP_lognormal_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rlnorm(1000*length(dose_g_4),log(means_g_4),sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("NP_lognormal_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rlnorm(1000*length(dose_e_4),log(means_e_4),sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("NP_lognormal_e_4_sim_%d.Rdata",i))
}
#inv-gaussian data
sd    <- c(rep(227.8176,2),rep(18528.14,2))
for (i in 1:4){
  means_g_5 <- X_g_5%*%t(iS[i,,drop=F]) + bkg[i]
  means_e_5 <- X_e_5%*%t(iS[i,,drop=F]) + bkg[i]
  means_g_4 <- X_g_4%*%t(iS[i,,drop=F]) + bkg[i]
  means_e_4 <- X_e_4%*%t(iS[i,,drop=F]) + bkg[i]
  
  doses = dose_g_5
  sim_data <- matrix(rinvgauss(1000*length(dose_g_5),means_g_5,sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("NP_invGaussian_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rinvgauss(1000*length(dose_e_5),means_e_5,sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("NP_invGaussian_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rinvgauss(1000*length(dose_g_4),means_g_4,sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("NP_invGaussian_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rinvgauss(1000*length(dose_e_4),means_e_4,sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("NP_invGaussian_e_4_sim_%d.Rdata",i))
}
setwd("..")