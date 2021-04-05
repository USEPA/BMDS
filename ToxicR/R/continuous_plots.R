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
cont_exp_3_f <-function(parms,d){
  g <- parms[1]
  b <- parms[2]
  e <- parms[4] 
  rval <- g*exp(-(b*d)^e)
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
    plot(dose,Response,type='n',...)
      
  }else{
    Response <- fit$data[,2,drop=F]
    doses = fit$data[,1,drop=F]
    plot(doses,Response,type='n',...)
  }
  

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
  temp_fit <- splinefun(test_doses,me)
  
  
  polygon(c(test_doses,test_doses[length(test_doses):1]),
          c(uq,lq[length(test_doses):1]),col = credint_col,border=credint_col)
  lines(test_doses,me)

  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    lines( c(fit$bmd[1],fit$bmd[1]),c(0,temp_fit(fit$bmd[1])))
    lines( c(fit$bmd[2],fit$bmd[2]),c(0,temp_fit(fit$bmd[2])))
    lines( c(fit$bmd[3],fit$bmd[3]),c(0,temp_fit(fit$bmd[3])))
  }
  
  if (BMD_DENSITY ==TRUE){
    temp = fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp = temp[!is.infinite(temp)]
    Dens =  density(temp,cut=c(max(test_doses)),adjust =1.5)
    Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.6
    temp = which(Dens$x < max(test_doses))
    D1_y = Dens$y[temp]
    D1_x = Dens$x[temp]
    qm = min(Response)
    polygon(c(0,D1_x,max(doses)),c(qm,qm+D1_y,qm),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
  }
  
  
  if (ncol(fit$data) ==4){
       points(doses,mean,...)
       arrows(x0=doses, y0=lerror, x1=doses, 
              y1=uerror, code=3, angle=90, length=0.1)
  }else{
       points(doses,Response,...)
  }
}
  

.plot.BMDcont_fit_maximized<-function(fit,qprob=0.05,...){
  
     density_col="blueviolet"
     credint_col="azure2"
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
     
    
  #sufficient statistics- This part dosen't makes senseArgument entry fixed 
  if (ncol(fit$data) == 4){ 
       mean <- fit$data[,2,drop=F]
       se   <- fit$data[,4,drop=F]/sqrt(fit$data[,3,drop=F])
       doses = fit$data[,1,drop=F]
       uerror <- mean+se
       lerror <- mean-se
       
       dose = c(doses,doses)
       Response = c(uerror,lerror)
       plot(dose,Response,type='n')
  }else{
    
       Response <- fit$data[,2,drop=F]
       doses = fit$data[,1,drop=F]
       plot(doses,Response,type='n')
  }
  # I fixed some logic of inputs in if/else statement- they used to be fit$data
  
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
  
  
  # This part should be loop 
  
  
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
  
  lines(test_doses,me)
  temp_fit <- splinefun(test_doses,me)
  
  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    lines( c(fit$bmd[1],fit$bmd[1]),c(0,temp_fit(fit$bmd[1])))
    lines( c(fit$bmd[2],fit$bmd[2]),c(0,temp_fit(fit$bmd[2])))
    lines( c(fit$bmd[3],fit$bmd[3]),c(0,temp_fit(fit$bmd[3])))
  }
  
  if (ncol(fit$data) == 4){
       points(doses,mean,...)
       arrows(x0=doses, y0=lerror, x1=doses, 
              y1=uerror, code=3, angle=90, length=0.1)
  }else{
       points(doses,Response,...)
  }
  
}

# Base plot- MCMC or BMD?
.plot.BMDcontinuous_MA <- function(A,qprob=0.05,...){
     density_col="blueviolet"
     credint_col="azure2"
     class_list <- names(A)
     
     fit_idx    <- grep("Individual_Model",class_list)
     
     #plot the model average curve
     if ("BMDcontinuous_MA_mcmc" %in% class(A)){ # mcmc run
          n_samps <- nrow(A[[fit_idx[1]]]$mcmc_result$PARM_samples); 
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
               plot(dose,Response,type='n')#,...)
               
          }else{
               Response <- data_d[,2,drop=F]
               doses = data_d[,1,drop=F]
               plot(jitter(doses),Response,type='n',...)
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
                    
                 
                    temp_f[ii,] <- cont_exp_3_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
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
          col1 = alphablend(credint_col,1)
          
          # Data structure for polygon - this part should be re-implmeneted as ggplot object
          polygon(c(test_doses,test_doses[length(test_doses):1]),
                  c(uq,lq[length(test_doses):1]),col = col1,border=col1)
         #test_dose = test_doses[is.finite(me)==TRUE]
         #me = me[is.finite(me) == TRUE]
          lines(test_doses,me,lwd=2)
          temp_fit <- splinefun(test_doses,me)
          bmd <- quantile(temp_bmd,c(qprob,0.5,1-qprob),na.rm = TRUE)
        
          lines( c(bmd[1],bmd[1]),c(0,temp_fit(bmd[1])))
          lines( c(bmd[2],bmd[2]),c(0,temp_fit(bmd[2])))
          lines( c(bmd[3],bmd[3]),c(0,temp_fit(bmd[3])))
          
          if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
               
               lines( c(bmd[1],bmd[1]),c(0,temp_fit(bmd[1])),lwd=2,lty=2)
               lines( c(bmd[2],bmd[2]),c(0,temp_fit(bmd[2])),lwd=3,)
               lines( c(bmd[3],bmd[3]),c(0,temp_fit(bmd[3])),lwd=2,lty=2)
          }
          
          
          temp = temp_bmd[!is.nan(temp_bmd)]
          temp = temp[!is.infinite(temp)]
          temp = temp[temp < 20 * max_dose]
          #return(temp)
          #print(c(max(temp),median(temp),min(temp)))
          Dens =  density(temp,cut=c(quantile(temp,0.995,na.rm = TRUE)),bw=10)
          Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.4
          temp = which(Dens$x < max(test_doses))
          D1_y = Dens$y[temp]
          D1_x = Dens$x[temp]
          qm = min(Response)
          polygon(c(0,D1_x,max(doses)),c(qm,qm+D1_y,qm),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
         
          for (ii in 1:length(fit_idx)){
               fit <- A[[fit_idx[ii]]]
               if (fit$model=="FUNL"){
                    f <- cont_FUNL_f(fit$fitted_model$parameters,test_doses)
               }  
               if (fit$model=="hill"){
                    f <- cont_hill_f(fit$fitted_model$parameters,test_doses)
               }
               if (fit$model=="exp-3"){
                   temp = fit$fitted_model$parameters 
                   temp = c(temp[1:2],0,temp[3],temp[4])
                    f <- cont_exp_3_f(temp,test_doses)
               }
               if (fit$model=="exp-5"){
                    f <- cont_exp_5_f(fit$fitted_model$parameters,test_doses)
               }
               if (fit$model=="power"){
                    f <- cont_power_f(fit$fitted_model$parameters,test_doses)
               }
               col = alphablend(col='coral3',min(1,A$posterior_probs[ii]*2))
               lines(test_doses,f,col=col,lwd = 2)
          }
          
     }else{
         
          data_d   <-  A[[fit_idx[1]]]$data
          max_dose <- max(data_d[,1])
          min_dose <- min(data_d[,1])
          test_doses <- seq(min_dose,max_dose,(max_dose-min_dose)/500); 
          temp_f   <- matrix(0,length(fit_idx),length(test_doses))
          
          if (ncol(data_d) == 4 ){ #sufficient statistics
               mean <- data_d[,2,drop=F]
               se   <- data_d[,4,drop=F]/sqrt(fit$data[,3,drop=F])
               doses = data_d[,1,drop=F]
               uerror <- mean+se
               lerror <- mean-se
               
               dose = c(doses,doses)
               Response = c(uerror,lerror)
               plot(dose,Response,type='n',...)
               
          }else{
               Response <- data_d[,2,drop=F]
               doses = data_d[,1,drop=F]
               plot(doses,Response,type='n',...)
          }
          
          for (ii in 1:length(fit_idx)){
               fit <- A[[fit_idx[ii]]]
               if (fit$model=="FUNL"){
                    temp_f[ii,] <- cont_FUNL_f(fit$parameters,test_doses)*A$posterior_probs[ii]
               }  
               if (fit$model=="hill"){
                    temp_f[ii,] <- cont_hill_f(fit$parameters,test_doses)*A$posterior_probs[ii]
               }
               if (fit$model=="exp-3"){
                    temp_f[ii,] <- cont_exp_3_f(fit$parameters,test_doses)*A$posterior_probs[ii]
               }
               if (fit$model=="exp-5"){
                    temp_f[ii,] <- cont_exp_5_f(fit$parameters,test_doses)*A$posterior_probs[ii]
               }
               if (fit$model=="power"){
                    temp_f[ii,] <- cont_power_f(fit$parameters,test_doses)*A$posterior_probs[ii]
               }
          }
          
          me <- colSums(temp_f)
          
          lines(test_doses,me,lwd=2)
         
          temp_fit <- splinefun(test_doses,me)
          bmds <- splinefun(A$ma_bmd[,2],A$ma_bmd[,1])
          temp_bmd <- bmds(runif(3000,0,max(A$ma_bmd[,2])))
          
          bmd <- quantile(temp_bmd,c(qprob,0.5,1-qprob),na.rm = TRUE)
          if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
               lines( c(bmd[1],bmd[1]),c(0,temp_fit(bmd[1])))
               lines( c(bmd[2],bmd[2]),c(0,temp_fit(bmd[2])))
               lines( c(bmd[3],bmd[3]),c(0,temp_fit(bmd[3])))
          }
          
          
          temp = temp_bmd[!is.nan(temp_bmd)]
          temp = temp[!is.infinite(temp)]
          Dens =  density(temp,cut=c(max(test_doses)),adjust =1.5)
          Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.4
          temp = which(Dens$x < max(test_doses))
          D1_y = Dens$y[temp]
          D1_x = Dens$x[temp]
          qm = min(Response)
          polygon(c(0,D1_x,max(doses)),c(qm,qm+D1_y,qm),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
          #plot the individual models proportional to their weight
          for (ii in 1:length(fit_idx)){
               fit <- A[[fit_idx[ii]]]
               if (fit$model=="FUNL"){
                    f <- cont_FUNL_f(fit$parameters,test_doses)
               }  
               if (fit$model=="hill"){
                    f <- cont_hill_f(fit$parameters,test_doses)
               }
               if (fit$model=="exp-3"){
                    f <- cont_exp_3_f(fit$parameters,test_doses)
               }
               if (fit$model=="exp-5"){
                    f <- cont_exp_5_f(fit$parameters,test_doses)
               }
               if (fit$model=="power"){
                    f <- cont_power_f(fit$parameters,test_doses)
               }
               col = alphablend(col='coral2',A$posterior_probs[ii])
               lines(test_doses,f,col=col,lwd = 2)
          }
          
          
     }
     
      
     if (ncol(fit$data) ==4){
          points(doses,mean)
          arrows(x0=doses, y0=lerror, x1=doses, 
                 y1=uerror, code=3, angle=90, length=0.1)
     }else{
          points(jitter(doses),Response,pch=16)
     }

}
