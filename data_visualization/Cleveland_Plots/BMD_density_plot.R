# Purpose: Plot stacked BMD density plots
# Author : Sooyeong Lim
# Date : 12/03/2020

# Load libraries 
library(ToxicR)
library(ggplot2)
library(tidyverse)
library(ggridges)

# Load required function from source code;
source("data_visualization/Cleveland_Plots/dicho_functions.R")

# Sample Dichotomous Data set
mData <- matrix(c(1, 3,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 40,50,
                  50, 48,50),nrow=5,ncol=3,byrow=T)

# Sample Fitting Case
A<-ma_dichotomous_fit(mData[,1],mData[,2],mData[,3], fit_type="mcmc", BMD_TYPE="added",BMR=0.1)






# Construct bmd sample plots for mcmc 
class_list <- names(A)
fit_idx    <- grep("Individual_Model",class_list)
qprob=0.05

#Dose levels
doses<-mData[,1]


for (i in fit_idx){
  # Loop for the model
  fit<-A[[i]]
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
  probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])
  doses = mData[,1,drop=T]
  
  if (fit$model=="hill"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_hill_f, d=test_doses)
    
  }
  if (fit$model=="gamma"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_gamma_f, d=test_doses)
    
  }
  if (fit$model=="logistic"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_logist_f, d=test_doses)
    
  }
  if (fit$model=="log-logistic"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_llogist_f, d=test_doses)
    
  }
  if (fit$model=="probit"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_probit_f, d=test_doses)
    
  }
  if (fit$model=="log-probit"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_lprobit_f, d=test_doses)
    
  }
  if (fit$model=="multistage"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_multistage_f, d=test_doses)
    
  }
  if (fit$model=="qlinear"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_qlinear_f, d=test_doses)
  }
  
  temp <- fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
  temp <- temp[!is.infinite(temp)]
  
  
  
  # Q <- t(Q)
  # 
  # me <- colMeans(Q)
  # lq <- apply(Q,2,quantile, probs = qprob)
  # uq <- apply(Q,2,quantile, probs = 1-qprob)
  
  Dens =  density(temp,cut=c(max(doses)))
  # what is this 0.4 means? Scale?
  
  # normalize it?-- We don't need it actually here
  # Dens$y = Dens$y/max(Dens$y) * max(probs)
  # temp = which(Dens$x < max(doses))
  # D1_y = Dens$y[temp]
  # D1_x = Dens$x[temp]
  
  
  
  # Do I need to stack up the dataset?
  
  
  temp_density<-data.frame(matrix(0,length(temp),3))
  temp_density[,2]=fit$model
  temp_density[,1]=temp
  temp_density[,3]=A$posterior_probs[i]
  
  assign(paste("t",i,sep="_"),temp_density)
  
}



# combine the fitting dataset here
t_combine<-rbind(t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9)


out_bmd<-ggplot()+
  geom_density(data=t_combine,aes(x=X1), fill = "blue", alpha=0.5)+
  labs(x="Dose Level (Dotted Line : MA BMD)",title="Cleveland Plot")+
  theme_classic()+xlim(c(min(doses),20))+
  facet_wrap(~X2,nrow=length(fit_idx), scales="free_y")+ geom_vline(xintercept = A$bmd[1],linetype="longdash")


out_bmd


samp<-t_combine %>% filter(X2=="hill")

q5<-quantile(samp$X1,0.05)
q95<-quantile(samp$X1,.95)
x.dens <- density(samp$X1)
df.dens <- data.frame(x = x.dens$x, y = x.dens$y)

# Object overwrapp
# I need to stack 8 different model. Is there any better way? free y_scale John's suggestion 
ggplot()+
  geom_density(data=samp,aes(x=X1),fill="blue",alpha=0.7)+
  xlim(c(0,20))+labs(x="Dose Level", title="BMD Estimates by Modles")+
  geom_area(data=subset(df.dens, x>=q5 &  x<=q95), aes(x=x,y=y),fill="blue",alpha=0.5)
# How can I add this for every case?


ggplot(data=t_combine,aes(x=X1, y=fct_reorder(X2,X3,.desc=T), fill = factor(stat(quantile))))+
  stat_density_ridges(geom="density_ridges_gradient",
                      calc_ecdf=TRUE,quantiles=c(0.025,0.975))+
  xlim(c(0,max(t_combine$X1)))+
  scale_fill_manual(name = "Probability", values = c("red", "blue", "red"), labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]"))+
  geom_vline(xintercept = A$bmd[1],linetype="longdash")+
  labs(y="Model",x="Dose Level (Dotted Line : MA BMD)", title="Density plots for each fitted model")+theme_classic()+
  theme(legend.position="none")






