#Create MA simulation
#Simulation Creation for the CMA paper
#Exp5.p copied from Simulation Conditions.Rmd (2/5/11)
library(ToxicR)
library(actuar)

set.seed(120302)
setwd('./Exp-5')
Exp5.p <- rbind(c(481,0.05,log(1/1.42870),2),
                c(481,0.02,log(1/1.42870),2),
                c(481,0.01,log(1/1.42870),2),
                c(481,0.1,log(1/1.42870),2) ,
                c(10.58,0.05,log(1.75),1.5),
                c(10.58,0.02,log(1.75),1.5),
                c(10.58,0.01,log(1.75),1.5),
                c(10.58,0.1, log(1.75),1.5))

dose_g_5 <- rep(c(0,6.25,12.5,25,50,100),each=10) 
dose_e_5 <- rep(c(0,20,40,60,80,100),each=10) 
dose_g_4 <- rep(c(0,12.5,25,50,100),each=10) 
dose_e_4 <- rep(c(0,25,50,75,100),each=10) 
sd    <- c(rep(77.5,4),rep(2.28,4))

doses2 <- seq(0,100,1)
lines(doses2,cont_exp_5_f(Exp5.p[2,],doses2))
#normal data
for (i in 1:8){
  means_g_5 <- cont_exp_5_f(Exp5.p[i,],dose_g_5)
  means_e_5 <- cont_exp_5_f(Exp5.p[i,],dose_e_5)
  means_g_4 <- cont_exp_5_f(Exp5.p[i,],dose_g_4)
  means_e_4 <- cont_exp_5_f(Exp5.p[i,],dose_e_4)
  
  doses = dose_g_5
  sim_data <- matrix(rnorm(1000*length(dose_g_5),means_g_5,sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_normal_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rnorm(1000*length(dose_e_5),means_e_5,sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_normal_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rnorm(1000*length(dose_g_4),means_g_4,sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_normal_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rnorm(1000*length(dose_e_4),means_e_4,sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_normal_e_4_sim_%d.Rdata",i))
}
#log-normal data
sd    <- c(rep(0.158,4),rep(0.209,4))
for (i in 1:8){
  means_g_5 <- cont_exp_5_f(Exp5.p[i,],dose_g_5)
  means_e_5 <- cont_exp_5_f(Exp5.p[i,],dose_e_5)
  means_g_4 <- cont_exp_5_f(Exp5.p[i,],dose_g_4)
  means_e_4 <- cont_exp_5_f(Exp5.p[i,],dose_e_4)
  
  doses = dose_g_5
  sim_data <- matrix(rlnorm(1000*length(dose_g_5),log(means_g_5),sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_lognormal_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rlnorm(1000*length(dose_e_5),log(means_e_5),sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_lognormal_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rlnorm(1000*length(dose_g_4),log(means_g_4),sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_lognormal_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rlnorm(1000*length(dose_e_4),log(means_e_4),sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_lognormal_e_4_sim_%d.Rdata",i))
}
#inv-gaussian data
sd    <- c(rep(18528.14,4),rep(227.8176,4))
#normal data
for (i in 1:8){
  means_g_5 <- cont_exp_5_f(Exp5.p[i,],dose_g_5)
  means_e_5 <- cont_exp_5_f(Exp5.p[i,],dose_e_5)
  means_g_4 <- cont_exp_5_f(Exp5.p[i,],dose_g_4)
  means_e_4 <- cont_exp_5_f(Exp5.p[i,],dose_e_4)
  
  doses = dose_g_5
  sim_data <- matrix(rinvgauss(1000*length(dose_g_5),means_g_5,sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_invGaussian_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rinvgauss(1000*length(dose_e_5),means_e_5,sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_invGaussian_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rinvgauss(1000*length(dose_g_4),means_g_4,sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_invGaussian_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rinvgauss(1000*length(dose_e_4),means_e_4,sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("Exp-5_invGaussian_e_4_sim_%d.Rdata",i))
}
setwd("..")