# Load ToxicR and ggplot package

library(tidyverse)
library(ggplot2)
library(ToxicR)
# Dichotomous case- Sample data input

mData <- matrix(c(1, 3,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 40,50,
                  50, 48,50),nrow=5,ncol=3,byrow=T)

# Sample RBMDS fitting
A<-ma_dichotomous_fit(mData[,1],mData[,2],mData[,3], fit_type="mcmc", BMR=0.1)

# Construct bmd sample plots for mcmc 
class_list <- names(A)

# Grap function extract # of indices from the text with same pattern
fit_idx    <- grep("Individual_Model",class_list)

# Create an empty matrix to contain BMD information from each model
bmd_ind<-matrix(0,length(fit_idx)+1,5)

for (i in 1:length(fit_idx)){
  
  # BMD -Median
  bmd_ind[i,1]<-A[[i]]$bmd[1]
  # BMD -5%
  bmd_ind[i,2]<-A[[i]]$bmd[2]
  # BMD -95%
  bmd_ind[i,3]<-A[[i]]$bmd[3]
  
  # Model name 
  bmd_ind[i,4]<-A[[i]]$model
  
  bmd_ind[i,5]<-A$posterior_probs[i]
  # 
  #   if (i==(length(fit_idx)+1)){
  #     
  #     # BMD -Median
  #     bmd_ind[i,1]<-A$bmd[1]
  #     # BMD -5%
  #     bmd_ind[i,2]<-A$bmd[2]
  #     # BMD -95%
  #     bmd_ind[i,3]<-A$bmd[3]
  #     
  #     # Model name 
  #     bmd_ind[i,4]<-"Model Average"
  #     
  #     bmd_ind[i,5]<-1
  #     
  #   }
  
  
}

bmd_ind[length(fit_idx)+1,1]<-A$bmd[1]
bmd_ind[length(fit_idx)+1,2]<-A$bmd[2]
bmd_ind[length(fit_idx)+1,3]<-A$bmd[3]

bmd_ind[length(fit_idx)+1,4]<-"Model Average"
bmd_ind[length(fit_idx)+1,5]<-1

bmd_ind_df<-data.frame(bmd_ind)

# The data structure will show 50% (BMD), 5% (BMDL), 95% upper bound from 9 different model
print(bmd_ind_df)

bmd_ind_df$X1
# I should add a vertical line from the bmd;

ggplot()+
  
  geom_point(data=bmd_ind_df, aes(x=as.numeric(X1), y=fct_reorder(X4,X5,.desc=T),size=(sqrt(as.numeric(X5)))), color="red")+
  #scale_colour_gradient(low = "gray", high = "black")+
  geom_vline(xintercept=as.numeric(bmd_ind_df$X1[10]), linetype="dotted")+
  theme_minimal()+
  labs(x="Dose Level",y="Model",title="BMD Estimates by Each Model (Sorted by Posterior Probability)",size="Posterior Probability")+
  theme(legend.position="none")+
  geom_errorbar(data=bmd_ind_df, width=0.2,aes(xmin=as.numeric(X2), xmax=as.numeric(X3), y=fct_reorder(X4,X5,.desc=T)),color="blue",alpha=0.3)






