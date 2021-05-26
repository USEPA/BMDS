# FUNL
cont_FUNL_f <- function(A,doses){
     b <- A[1] + A[2]*exp(-exp(A[6])*(doses-A[5])^2)*(1/(1+exp(-(doses-A[3])/A[4])))
     return(b)
}

#dichotomous hill
cont_hill_f <- function(parms,d){
  g  <- parms[1] 
  nu <- parms[2]
  k  <- parms[3];
  n  <- parms[4]; 
  rval <- g + nu*d^n/(k^n+d^n)
  return (rval)
}
#dichotomous log-logistic
cont_exp_5_f <- function(parms,d){
  g <- parms[1]
  b <- parms[2];
  c <- parms[3];
  e <- parms[4]; 
  rval <- g*(exp(c)-(exp(c)-1.0)*(exp(-(b*d)^e)))
  return (rval)
}

#
cont_exp_3_f <-function(parms,d,decrease = TRUE){
  if (decrease){
    f_sign = -1; 
  }else{
    f_sign = 1; 
  }
  g <- parms[1]
  b <- parms[2]
  e <- parms[3] 
  rval <- g*exp(f_sign*(b*d)^e)
  return (rval)
}

cont_power_f <-function(parms,d){
  g <- parms[1]; 
  b <- parms[2];
  a <- parms[3]; 
  rval <- g + b*d^a
  return (rval)
}



.plot.BMDcont_fit_MCMC<-function(fit,qprob=0.05,...){
  
  fit<-A
  density_col="blueviolet"
  credint_col="azure2"
  BMD_DENSITY = T
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
  
  
  
  if (ncol(fit$data) == 4 ){ #sufficient statistics
    mean <- fit$data[,2,drop=F]
    se   <- fit$data[,4,drop=F]/sqrt(fit$data[,3,drop=F])
    doses = fit$data[,1,drop=F]
    uerror <- mean+se
    lerror <- mean-se
    
    dose = c(doses,doses)
    Response = c(uerror,lerror)
    # plot(dose,Response,type='n',...)
  }else{
    Response <- fit$data[,2,drop=F]
    doses = fit$data[,1,drop=F]
    # plot(doses,Response,type='n',...)
  }
  
  # Single Model 
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
  
  if (fit$model=="FUNL"){
     Q <- apply(fit$mcmc_result$PARM_samples,1,cont_FUNL_f, d=test_doses)   
  }
  if (fit$model=="hill"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,cont_hill_f, d=test_doses)
  }
  if (fit$model=="exp-3"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,cont_exp_3_f, d=test_doses)
  }
  if (fit$model=="exp-5"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,cont_exp_5_f, d=test_doses)
  }
  if (fit$model=="power"){
    Q <- apply(fit$mcmc_result$PARM_samples,1,cont_power_f, d=test_doses)
  }
  
 
  Q <- t(Q)
  me <- colMeans(Q)
  lq <- apply(Q,2,quantile, probs = qprob)
  uq <- apply(Q,2,quantile, probs = 1-qprob)
  
  # Continous case density? 
  temp_fit <- splinefun(test_doses,me)
  
  
  
  # Geom_polygon ? etc..
  
  out<-ggplot()+
    geom_line(aes(x=test_doses,y=me),color="blue")+
    labs(x="Dose", y="Response",title=paste(fit$fitted_model$full_model, "MCMC",sep=",  Fit Type: " ))+theme_minimal()

  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    out2<-out+
      geom_segment(aes(x=fit$bmd, y=temp_fit(x=fit$bmd), xend=fit$bmd, yend=0),color="Red")
  }
  
# Add density 
  if (BMD_DENSITY ==TRUE){
    temp = fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp = temp[!is.infinite(temp)]
    Dens =  density(temp,cut=c(max(test_doses)),adjust =1.5)
    Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.6
    temp = which(Dens$x < max(test_doses))
    D1_y = Dens$y[temp]
    D1_x = Dens$x[temp]
    qm = min(Response)
    
    
    out3<-out2+geom_polygon(aes(x=c(0,D1_x,max(doses)),y=c(0,0+D1_y,0)), fill = "lightblue1", alpha=0.6)
    # polygon(c(0,D1_x,max(doses)),c(qm,qm+D1_y,qm),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
  }
  
  # 
  # if (ncol(fit$data) ==4){
  #      points(doses,mean,...)
  #      arrows(x0=doses, y0=lerror, x1=doses, 
  #             y1=uerror, code=3, angle=90, length=0.1)
  # }else{
  #      points(doses,Response,...)
  # }
  # 
  
  out4<-out3+geom_point(aes(x=doses,y=Response))

  out5<-out4+geom_polygon(aes(x=c(test_doses,test_doses[length(test_doses):1]),y=c(uq,lq[length(test_doses):1])),fill="blue",alpha=0.1)
  out5
  
}
  

.plot.BMDcont_fit_maximized<-function(A,qprob=0.05,...){
  
  
  fit<-A
  density_col="blueviolet"
  credint_col="azure2"
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
     
  
  if (ncol(fit$data) == 4 ){ #sufficient statistics
    mean <- fit$data[,2,drop=F]
    se   <- fit$data[,4,drop=F]/sqrt(fit$data[,3,drop=F])
    doses = fit$data[,1,drop=F]
    uerror <- mean+se
    lerror <- mean-se
    dose = c(doses,doses)
    Response = c(uerror,lerror)
    # plot(dose,Response,type='n',...)
  }else{
    Response <- fit$data[,2,drop=F]
    doses = fit$data[,1,drop=F]
    # plot(doses,Response,type='n',...)
  }
 
  
  # I fixed some logic of inputs in if/else statement- they used to be fit$data
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)-min(doses))/300)
  
  
  if (fit$model=="FUNL"){
     me <- cont_FUNL_f(fit$parameters,test_doses)
  }  
  if (fit$model=="hill"){
    me <- cont_hill_f(fit$parameters,test_doses)
  }
  if (fit$model=="exp-3"){
    me <- cont_exp_3_f(fit$parameters,test_doses)
  }
  if (fit$model=="exp-5"){
    me <- cont_exp_5_f(fit$parameters,test_doses)
  }
  if (fit$model=="power"){
    me <- cont_power_f(fit$parameters,test_doses)
  }
  

  temp_fit <- splinefun(test_doses,me)
  
  
  out<-ggplot()+
          geom_line(aes(x=test_doses,y=me),color="blue")+
          labs(x="Dose", y="Response",title=paste(fit$full_model, "Maximized",sep=",  Fit Type: " ))+
          theme_minimal() + ylim(c(min(Response)*0.95,max(Response)*1.05)) +
          xlim(c(min(test_doses),max(test_doses)))
                   
  
  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    out2<-out+
      geom_segment(aes(x=fit$bmd, y=temp_fit(x=fit$bmd), xend=fit$bmd, yend=min(Response)*0.95),color="Red")
  }
  
  data_in<-data.frame(cbind(doses,Response))
  out3<-out2+
        geom_point(data=data_in,aes(x=Dose,y=Resp))
  
  
  out3
  
}


# Base plot- MCMC or BMD?
.plot.BMDcontinuous_MA <- function(A,qprob=0.05,...){
  
  # Should be matched with BMD_MA plots
    
     density_col="blueviolet"
     credint_col="azure2"
     class_list <- names(A)
     
     fit_idx    <- grep("Individual_Model",class_list)
     
     #plot the model average curve
     if ("BMDcontinuous_MA_mcmc" %in% class(A)){ # mcmc run
          n_samps <- nrow(A[[fit_idx[1]]]$mcmc_result$PARM_samples)
          data_d   <-  A[[fit_idx[1]]]$data
          max_dose <- max(data_d[,1])
          min_dose <- min(data_d[,1])
          test_doses <- seq(min_dose,max_dose,(max_dose-min_dose)/200); 
          ma_samps <- sample(fit_idx,n_samps, replace=TRUE,prob = A$posterior_probs)
          temp_f   <- matrix(0,n_samps,length(test_doses))
          temp_bmd <- rep(0,length(test_doses))
          
          
          if (ncol(data_d) == 4 ){ #sufficient statistics
            mean <- data_d[,2,drop=F]
            se   <- data_d[,4,drop=F]/sqrt(fit$data[,3,drop=F])
            doses = data_d[,1,drop=F]
            uerror <- mean+se
            lerror <- mean-se
            dose = c(doses,doses)
            Response = c(uerror,lerror)
            lm_fit = lm(mean = doses,weights = 1/se*se)
          }else{
            Response <- data_d[,2,drop=F]
            doses = data_d[,1,drop=F]
            lm_fit = lm(Response~doses)
          }
          
          if (coefficients(lm_fit)[2] < 0){
            decrease = TRUE
          }else{
            decrease = FALSE
          }
          
          for (ii in 1:n_samps){
               fit <- A[[fit_idx[ma_samps[ii]]]]
               if (fit$model=="FUNL"){
                    temp_f[ii,] <- cont_FUNL_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }  
               if (fit$model=="hill"){
                    temp_f[ii,] <- cont_hill_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
               if (fit$model=="exp-3"){

                    temp_f[ii,] <- cont_exp_3_f(fit$mcmc_result$PARM_samples[ii,],test_doses,decrease)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
               if (fit$model=="exp-5"){
                    temp_f[ii,] <- cont_exp_5_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
               if (fit$model=="power"){
                    temp_f[ii,] <- cont_power_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
                    temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
               }
          }
          temp_f[is.infinite(temp_f)] = NA
        
          me <- colMeans(temp_f,na.rm = TRUE)
          lq <- apply(temp_f,2,quantile, probs = qprob,na.rm = TRUE)
          uq <- apply(temp_f,2,quantile, probs = 1-qprob,na.rm = TRUE)
 
          
         
          plot_gg<-ggplot()+
                  geom_point(aes(x=doses,y=Response))+
                  xlim(c(min(doses),max(doses)*1.03))+
                  labs(x="Dose", y="Proportion",title="Continous MA fitting")+
                  theme_minimal()
          
         plot_gg<-plot_gg+
                  geom_ribbon(aes(x=test_doses,ymin=lq,ymax=uq),fill="blue",alpha=0.1)
         
         plot_gg<-plot_gg+
                  geom_line(aes(x=test_doses,y=me),col="blue",size=2)
         
          bmd <- quantile(temp_bmd,c(qprob,0.5,1-qprob),na.rm = TRUE)

          ## Plot the CDF of the Posterior
          if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
            temp = temp_bmd[temp_bmd < 10*max(test_doses)]
            temp = temp[!is.infinite(temp_bmd)]
            temp = temp[!is.na(temp)]
       
            Dens =  density(temp,cut=c(max(test_doses)))
          
            Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.6
            temp = which(Dens$x < max(test_doses))
            D1_y = Dens$y[temp]
            D1_x = Dens$x[temp]
            qm = min(Response)
            scale = (max(Response)-min(Response))/max(D1_y) *.75
            # BMD MA density needs to be double checked 
            plot_gg<-plot_gg+
                    geom_polygon(aes(x=c(max(0,min(D1_x)),D1_x,max(0,min(D1_x))),
                                     y=c(min(Response),min(Response)+D1_y*scale,min(Response))),
                                     fill = "blueviolet", alpha=0.6)

           }
          ## 
          # Add lines to the BMD
          ma_mean <- splinefun(test_doses,me)
          ma_BMD = A$bmd
          plot_gg = plot_gg + 
                     geom_segment(aes(x=A$bmd, y=ma_mean(A$bmd), xend=A$bmd, yend=min(Response)),color="Red")
          
           
          #Plot only level >2

          df<-NULL
         
          for (ii in 1:length(fit_idx)){
            
            if (A$posterior_probs[ii]>0.05){
               fit <- A[[fit_idx[ii]]]
               if (fit$model=="FUNL"){
                    f <- cont_FUNL_f(fit$fitted_model$parameters,test_doses)
               }  
               if (fit$model=="hill"){
                    f <- cont_hill_f(fit$fitted_model$parameters,test_doses)
               }
               if (fit$model=="exp-3"){
                   temp = fit$fitted_model$parameters
                    f <- cont_exp_3_f(temp,test_doses,decrease)
               }
               if (fit$model=="exp-5"){
                    f <- cont_exp_5_f(fit$fitted_model$parameters,test_doses)
               }
               if (fit$model=="power"){
                    f <- cont_power_f(fit$fitted_model$parameters,test_doses)
               }
               col = alphablend(col='coral3',A$posterior_probs[ii])
               # Not using loop, but save data in the external data and load it later
               temp_df<-data.frame(x_axis=test_doses,y_axis=f,cols=col,model_no=ii, alpha_lev=A$posterior_probs[ii])
               df<-rbind(df,temp_df)
               # Not using loop, but save data in the external data and load it later
               if (A$posterior_probs[ii] > 0.01){
                   temp_df<-data.frame(x_axis=test_doses,y_axis=f,cols=col,model_no=ii, alpha_lev=A$posterior_probs[ii])
                   df<-rbind(df,temp_df)
               }
            }
            
          }
          
          plot_gg<- plot_gg+
                 geom_line(data=df, aes(x=x_axis,y=y_axis,color=cols),alpha=0.5,show.legend=F)+
                 theme_minimal()

     }
     else{ # mcmc run
       
       data_d   <-  A[[fit_idx[1]]]$data
       max_dose <- max(data_d[,1])
       min_dose <- min(data_d[,1])
       test_doses <- seq(min_dose,max_dose,(max_dose-min_dose)/200); 
       temp_bmd <- rep(0,length(test_doses))
       
       if (ncol(data_d) == 4 ){ #sufficient statistics
         mean <- data_d[,2,drop=F]
         se   <- data_d[,4,drop=F]/sqrt(fit$data[,3,drop=F])
         doses = data_d[,1,drop=F]
         uerror <- mean+se
         lerror <- mean-se
         dose = c(doses,doses)
         Response = c(uerror,lerror)
         lm_fit = lm(mean = doses,weights = 1/se*se)
       }else{
         Response <- data_d[,2,drop=F]
         doses = data_d[,1,drop=F]
         lm_fit = lm(Response~doses)
       }
       
       
       if (coefficients(lm_fit)[2] < 0){
         decrease = TRUE
       }else{
         decrease = FALSE
       }
       me = test_doses*0   
       for (ii in 1:length(fit_idx)){
         fit <- A[[fit_idx[ii]]]
         if (fit$model=="FUNL"){
           t <- cont_FUNL_f(fit$parameters,test_doses)
           if(BB$posterior_probs[ii] > 0){
             me = t*BB$posterior_probs[ii] + me
           }
          
         }  
         if (fit$model=="hill"){
            
           t <- cont_hill_f(fit$parameters,test_doses)
           if(BB$posterior_probs[ii] > 0){
             me = t*BB$posterior_probs[ii] + me
           }
         }
         if (fit$model=="exp-3"){
           t <- cont_exp_3_f(fit$parameters,test_doses,decrease)
   
           if(BB$posterior_probs[ii] > 0){
             me = t*BB$posterior_probs[ii] + me
           }
         }
         if (fit$model=="exp-5"){
           t <- cont_exp_5_f(fit$parameters,test_doses)
           if(BB$posterior_probs[ii] > 0){
             me = t*BB$posterior_probs[ii] + me
           }
         }
         if (fit$model=="power"){
           t <- cont_power_f(fit$parameters,test_doses)
           if(BB$posterior_probs[ii] > 0){
             me = t*BB$posterior_probs[ii] + me
           }
         }
       }

       plot_gg<-ggplot()+
         geom_point(aes(x=doses,y=Response))+
         xlim(c(min(doses),max(doses)*1.03))+
         labs(x="Dose", y="Proportion",title="Continous MA fitting")+
         theme_minimal()
       
        
       plot_gg<-plot_gg+
         geom_line(aes(x=test_doses,y=me),col="blue",size=2)
       
      
       ## 
       # Add lines to the BMD
       ma_mean <- splinefun(test_doses,me)
       ma_BMD = A$bmd
       plot_gg = plot_gg + 
         geom_segment(aes(x=A$bmd, y=ma_mean(A$bmd), xend=A$bmd, yend=min(Response)),color="Red")
       
       
       #Plot only level >2
       
       df<-NULL
       for (ii in 1:length(fit_idx)){
         
         if (A$posterior_probs[ii]>0.05){
           fit <- A[[fit_idx[ii]]]
           if (fit$model=="FUNL"){
             f <- cont_FUNL_f(fit$parameters,test_doses)
           }  
           if (fit$model=="hill"){
             f <- cont_hill_f(fit$parameters,test_doses)
           }
           if (fit$model=="exp-3"){
             temp = fit$parameters 
             #temp = c(temp[1:2],0,temp[3],temp[4])
             f <- cont_exp_3_f(temp,test_doses,decrease)
           }
           if (fit$model=="exp-5"){
             f <- cont_exp_5_f(fit$parameters,test_doses)
           }
           if (fit$model=="power"){
             f <- cont_power_f(fit$parameters,test_doses)
           }
           
           col = alphablend(col='coral3',A$posterior_probs[ii])
           # Not using loop, but save data in the external data and load it later
           temp_df<-data.frame(x_axis=test_doses,y_axis=f,cols=col,model_no=ii, alpha_lev=A$posterior_probs[ii])
           df<-rbind(df,temp_df)
         }
       }
       
       plot_gg<- plot_gg+
         geom_line(data=df, aes(x=x_axis,y=y_axis,color=cols),alpha=0.5,show.legend=F)+
         theme_minimal()
       
     }
     return(plot_gg)
}
