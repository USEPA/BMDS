#Create MA simulation
#Simulation Creation for the CMA paper
library(ToxicR)
library(actuar)

set.seed(101711)
setwd('./Hill')
Hill.p <- rbind(c(481,-144.3,70,3.3),
                c(481,-144.3,40,1.3),
                c(481,-144.2,15,1.1),
                c(481,-144.3,50,4) ,
                c(10.58,9.7,70,3.5),
                c(10.58,9.7,25,3),
                c(10.58,9.7,15,2),
                c(10.58,9.7,50,4))

dose_g_5 <- rep(c(0,6.25,12.5,25,50,100),each=10) 
dose_e_5 <- rep(c(0,20,40,60,80,100),each=10) 
dose_g_4 <- rep(c(0,12.5,25,50,100),each=10) 
dose_e_4 <- rep(c(0,25,50,75,100),each=10) 
sd    <- c(rep(77.5,4),rep(2.28,4))
#normal data
for (i in 1:8){
  means_g_5 <- cont_hill_f(Hill.p[i,],dose_g_5)
  means_e_5 <- cont_hill_f(Hill.p[i,],dose_e_5)
  means_g_4 <- cont_hill_f(Hill.p[i,],dose_g_4)
  means_e_4 <- cont_hill_f(Hill.p[i,],dose_e_4)
  
  doses = dose_g_5
  sim_data <- matrix(rnorm(1000*length(dose_g_5),means_g_5,sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_normal_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rnorm(1000*length(dose_e_5),means_e_5,sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_normal_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rnorm(1000*length(dose_g_4),means_g_4,sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_normal_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rnorm(1000*length(dose_e_4),means_e_4,sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_normal_e_4_sim_%d.Rdata",i))
}
#log-normal data
sd    <- c(rep(0.158,4),rep(0.209,4))
for (i in 1:8){
  means_g_5 <- cont_hill_f(Hill.p[i,],dose_g_5)
  means_e_5 <- cont_hill_f(Hill.p[i,],dose_e_5)
  means_g_4 <- cont_hill_f(Hill.p[i,],dose_g_4)
  means_e_4 <- cont_hill_f(Hill.p[i,],dose_e_4)
  
  doses = dose_g_5
  sim_data <- matrix(rlnorm(1000*length(dose_g_5),log(means_g_5),sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_lognormal_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rlnorm(1000*length(dose_e_5),log(means_e_5),sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_lognormal_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rlnorm(1000*length(dose_g_4),log(means_g_4),sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_lognormal_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rlnorm(1000*length(dose_e_4),log(means_e_4),sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_lognormal_e_4_sim_%d.Rdata",i))
}
#inv-gaussian data
sd    <- c(rep(18528.14,4),rep(227.8176,4))
#normal data
for (i in 1:8){
  means_g_5 <- cont_hill_f(Hill.p[i,],dose_g_5)
  means_e_5 <- cont_hill_f(Hill.p[i,],dose_e_5)
  means_g_4 <- cont_hill_f(Hill.p[i,],dose_g_4)
  means_e_4 <- cont_hill_f(Hill.p[i,],dose_e_4)
  
  doses = dose_g_5
  sim_data <- matrix(rinvgauss(1000*length(dose_g_5),means_g_5,sd[i]),nrow=1000,ncol=length(dose_g_5),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_invGaussian_g_5_sim_%d.Rdata",i))
  doses = dose_e_5
  sim_data <- matrix(rinvgauss(1000*length(dose_e_5),means_e_5,sd[i]),nrow=1000,ncol=length(dose_e_5),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_invGaussian_e_5_sim_%d.Rdata",i))
  doses = dose_g_4
  sim_data <- matrix(rinvgauss(1000*length(dose_g_4),means_g_4,sd[i]),nrow=1000,ncol=length(dose_g_4),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_invGaussian_g_4_sim_%d.Rdata",i))
  doses = dose_e_4
  sim_data <- matrix(rinvgauss(1000*length(dose_e_4),means_e_4,sd[i]),nrow=1000,ncol=length(dose_e_4),byrow=T)
  save(doses,sim_data,file=sprintf("Hill_invGaussian_e_4_sim_%d.Rdata",i))
}