#Set the default clevland_plot method generic for all of the classes. 
MAdensity_plot <- function (A, ...){
  source("dicho_functions.R")
  UseMethod("MAdensity_plot")
}

# Sample Dichotomous Data set


.plot.density.BMDdichotomous_MA_MCMC<-function(A){
# Construct bmd sample plots for mcmc
  class_list <- names(A)
  fit_idx    <- grep("Individual_Model",class_list)
  qprob=0.05
  
  #Dose levels
  data<-A$Individual_Model_1$data
  doses<-data[,1]
  
  
  for (i in fit_idx){
    # Loop for the model
    fit<-A[[i]]
    test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
    probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])

    
    
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
    
  # Combine the fitting dataset here
  t_combine<-rbind(t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9)
  
  #Plot 
  ggplot(data=t_combine,aes(x=X1, y=fct_reorder(X2,X3,.desc=T), fill = factor(stat(quantile))))+
    stat_density_ridges(geom="density_ridges_gradient",calc_ecdf=TRUE,quantiles=c(0.025,0.975))+
    xlim(c(0,max(t_combine$X1)))+
    scale_fill_manual(name = "Probability", values = c("red", "blue", "red"), 
                      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]"))+
    geom_vline(xintercept = A$bmd[1],linetype="longdash")+
    labs(y="Model",x="Dose Level (Dotted Line : MA BMD)", title="Density plots for each fitted model (Fit type: MCMC)")+theme_classic()+
    theme(legend.position="none")

}

# Define another function for mle/laplace
.plot.density.BMDdichotomous_MA_maximized<-function(A,fit_type="laplace"){
  
  
  
  class_list <- names(A)
  
  if (fit_type=="laplace"){
    fit_idx <- grep("Fitted_Model",class_list)
    qprob=0.05
    
    #Dose levels
    data<-A$Fitted_Model_1$data
    doses<-data[,1]
    
  }
  
  else{
  fit_idx <- grep("Individual_Model",class_list)
  qprob=0.05
  
  #Dose levels
  data<-A$Individual_Model_1$data
  doses<-data[,1]
  }
  
  for (i in fit_idx){
    fit<-A[[i]]
    test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
    
    if (fit$model=="hill"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      #fit$parameters[2] = .logit(fit$parameters[2])
      me <- .dich_hill_f(fit$parameters, d=test_doses)
    }
    if (fit$model=="gamma"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      me <- .dich_gamma_f(fit$parameters, d=test_doses)
    }
    if (fit$model=="logistic"){
      me <- .dich_logist_f(fit$parameters, d=test_doses)
    }
    if (fit$model=="log-logistic"){
      #fit$parameters[1] = logit(fit$parameters[1])
      me <- .dich_llogist_f(fit$parameters, d=test_doses)    
    }
    if (fit$model=="probit"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      print
      me <- .dich_probit_f(fit$parameters, d=test_doses)    
    }
    if (fit$model=="log-probit"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      me <- .dich_lprobit_f(fit$parameters, d=test_doses)    
    }
    if (fit$model=="multistage"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      me <- .dich_multistage_f(fit$parameters, d=test_doses)    
    }
    if (fit$model=="qlinear"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      me <- .dich_qlinear_f(fit$parameters, d=test_doses)    
    }
    if (fit$model=="weibull"){
      #fit$parameters[1] = .logit(fit$parameters[1])
      me <- .dich_weibull_f(fit$parameters, d=test_doses)    
    }

    temp <- fit$bmd_dist[,1]
    temp <- temp[!is.infinite(temp)]
    
    
    
    temp_density<-data.frame(matrix(0,length(temp),3))
    temp_density[,2]=fit$model
    temp_density[,1]=temp
    temp_density[,3]=A$posterior_probs[i]
    
    assign(paste("t",i,sep="_"),temp_density)
    
  }
  
  t_combine<-rbind(t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9)
  
  #Plot 
  ggplot(data=t_combine,aes(x=X1, y=fct_reorder(X2,X3,.desc=T), fill = factor(stat(quantile))))+
  stat_density_ridges(geom="density_ridges_gradient",calc_ecdf=TRUE,quantiles=c(0.025,0.975))+
  xlim(c(0,max(t_combine$X1)))+
    scale_fill_manual(name = "Probability", values = c("red", "blue", "red"), 
                      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]"))+
    geom_vline(xintercept = quantile(A$BMD_CDF[,1],0.5),linetype="longdash")+
    labs(y="Model",x="Dose Level (Dotted Line : MA BMD)", title="Density plots for each fitted model (Fit type: Maximized)")+theme_classic()+
    theme(legend.position="none")
  
  
}




