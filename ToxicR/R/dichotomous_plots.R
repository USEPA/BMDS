# Dichotomous functions are defined here
{
  .logit <- function(p)
  {
    return (log(p/(1-p)))
  }
  
  #dichotomous hill
  .dich_hill_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    n <- 1/(1+exp(-parms[2])); 
    a <- parms[3];
    b <- parms[4]; 
    rval <- g + (1-g)*n*(1/(1+exp(-a-b*log(d))))
    return (rval)
  }
  #dichotomous log-logistic
  .dich_llogist_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
    return (rval)
  }
  #dichotomous log-probit
  .dich_lprobit_f <-function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
    return (rval)
  }
  
  #dichotomous weibull
  .dich_weibull_f <-function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*(1-exp(-b*d^a))
    return (rval)
  }
  
  #dichotomous gamma
  .dich_gamma_f <-function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    b <- parms[3]; 
    rval <- g + (1-g)*pgamma(b*d,a,1)
    return (rval)
  }
  
  #dichtomous logistic
  .dich_logist_f <- function(parms,d){
    rval <- 1/(1+exp(-parms[1]-parms[2]*d))
    return (rval)
  }
  
  #dichtomous probit
  .dich_probit_f <- function(parms,d){
    rval <- pnorm(parms[1]+parms[2]*d)
    return (rval)
  }
  
  .dich_qlinear_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    a <- parms[2];
    return (g + (1-g)*1-exp(-a*d))
  }
  
  .dich_multistage_f <- function(parms,d){
    g <- 1/(1+exp(-parms[1])); 
    rval = d*0
    for (ii  in 2:length(parms)){
      rval = rval - parms[ii]*d^(ii-1)
    }
    return (g + (1-g)*1-exp(rval))
  }
  
  
}

{
  .plot.BMDdich_fit_MCMC <-function(fit,fit_type,qprob=0.05,...){
    
    density_col="red"
    credint_col="lightblue1"
    BMD_DENSITY = T
    
    if (qprob < 0 || qprob > 0.5){
      stop( "Quantile probability must be between 0 and 0.5")
    }
    
    
    
    # How this is calculated?
    # This part - how it is derived?
    probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])
    se <- sqrt(probs*(1-probs)/fit$data[,3,drop=T])
    
    
    doses = fit$data[,1,drop=T]
    uerror <- apply(cbind(probs*0+1,probs+se),1,min)
    lerror <- apply(cbind(probs*0,probs-se),1,max)
    
    dose = c(doses,doses)
    Response = c(uerror,lerror)
    
    # Basic structure of display
    # main should show the models' information, I think this aprt sohuld be fixed.
    
    # Dichotomous's response is between 0 to 1
    
    # Change this to ggplot object
    
    
    #plot(dose,Response,type='n',main=fit$fitted_model$full_model)
    
    # We need to adjust the range here too
    
    # S3 object not fitted here for the title part
    out<-ggplot()+
      geom_errorbar(aes(x=doses, ymin=lerror, ymax=uerror),color="grey")+xlim(c(min(dose),max(dose)))+ylim(c(0,1))+labs(x="Dose", y="Proportion",title=paste(fit$fitted_model$full_model, fit_type,sep=",  Fit Type: " ))+theme_minimal()
    
    
    
    
    test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/100)
    
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
    
    if (fit$model=="weibull"){
      Q <- apply(fit$mcmc_result$PARM_samples,1,.dich_weibull_f, d=test_doses)
      
    }
    
    
    temp <- fit$mcmc_result$BMD_samples[!is.nan(fit$mcmc_result$BMD_samples)]
    temp <- temp[!is.infinite(temp)]
    test <- density(temp)
    
    
    Q <- t(Q)
    
    me <- colMeans(Q)
    lq <- apply(Q,2,quantile, probs = qprob)
    uq <- apply(Q,2,quantile, probs = 1-qprob)
    
    # Splien function is used to test column everage from MCMC
    
    temp_fit <- splinefun(test_doses,me)
    
    
    # Object 2
    # Option- color , border- credint_col 
    # Ploygon object alpha need to be decided- should be transparent
    
    # Input - color , check 
    #polygon(c(test_doses,test_doses[length(test_doses):1]),
    #        c(uq,lq[length(test_doses):1]),col = credint_col, border=credint_col)
    
    # Polygon changed to Geom_ribbon
    out2<-out+geom_ribbon(aes(x=test_doses,ymin=lq,ymax=uq),fill="blue",alpha=0.1)
    out2<-out+geom_polygon(aes(x=c(test_doses[length(test_doses):1],test_doses),y=c(uq[length(test_doses):1],lq)),fill="blue",alpha=0.1)
    
    
    out3<-out2+geom_smooth(aes(x=test_doses,y=me),col="blue")+geom_point(aes(x=doses,y=probs))
    
    
    # This part is for referecne 
    #geom_segment(data=bmd_dots, aes(x=H.bmd, y=.dich_weibull_f.H.fitted_model.parameters..H.bmd., xend=H.bmd, yend=0), color="Red")+
    
    out4<-out3+geom_segment(aes(x=fit$bmd, y=temp_fit(x=fit$bmd), xend=fit$bmd, yend=0), color="Red")
    out4
    
    
    
    
    # 
    # if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    #   lines( c(fit$bmd[1],fit$bmd[1]),c(0,temp_fit(fit$bmd[1])))
    #   lines( c(fit$bmd[2],fit$bmd[2]),c(0,temp_fit(fit$bmd[2])))
    #   lines( c(fit$bmd[3],fit$bmd[3]),c(0,temp_fit(fit$bmd[3])))
    # }
    # 
    
    # Adding density 
    
    # Object 3 - Density object - In Shiny we can on / off this 
    # Density - Polygon/Other option?
    if (BMD_DENSITY ==TRUE){
      Dens =  density(temp,cut=c(max(doses)))
      # what is this 0.4 means? Scale?
      Dens$y = Dens$y/max(Dens$y) * max(probs)*0.4
      temp = which(Dens$x < max(doses))
      D1_y = Dens$y[temp]
      D1_x = Dens$x[temp]
      
      #geom ploygon 
      # geom_polygon(aes(x=c(test_doses,test_doses[length(test_doses):1]),y=c(uq,lq[length(test_doses):1])), fill="blue",alpha=0.1)
      #polygon(c(0,D1_x,max(doses)),c(0,D1_y,0),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
      
      out5<-out4+geom_polygon(aes(x=c(0,D1_x,max(doses)),y=c(0,D1_y,0)), fill = "lightblue1", alpha=0.5)
      
      out5
    }
    
    
    out6<-ggplotly(out5)
    # Already reflected from above 
    #points(doses,probs)
    #arrows(x0=doses, y0=lerror, x1=doses, 
    #       y1=uerror, code=3, angle=90, length=0.1)
    
    out6
  }
  



.plot.BMDdich_fit_maximized <- function(fit,qprob=0.05,...){
  
  density_col="red"
  credint_col="lightblue1"
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
   
  
  probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])
  se <- sqrt(probs*(1-probs)/fit$data[,3,drop=T])
  doses = fit$data[,1,drop=T]
  uerror <- apply(cbind(probs*0+1,probs+2*se),1,min)
  lerror <- apply(cbind(probs*0,probs-2*se),1,max)
  
  dose = c(doses,doses)
  Response = c(uerror,lerror)
  plot(dose,Response,type='n',main=fit$full_model...)
  
  
  test_doses <- seq(min(doses),max(doses)*1.03,(max(doses)*1.03-min(doses))/300)
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
  
 
  lines(test_doses,me)
  temp_fit <- splinefun(test_doses,me)
  
  if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
    lines( c(fit$bmd[1],fit$bmd[1]),c(0,temp_fit(fit$bmd[1])))
    lines( c(fit$bmd[2],fit$bmd[2]),c(0,temp_fit(fit$bmd[2])))
    lines( c(fit$bmd[3],fit$bmd[3]),c(0,temp_fit(fit$bmd[3])))
  }
  temp <- fit$bmd_dist
  temp <- temp[!is.infinite(temp[,1]),]
  temp <- temp[!is.nan(temp[,1]),]
  cdf_rows = nrow(temp)
  a = temp[,1]
  b = temp[,2]
  b = b[(!is.infinite(a))*(!is.na(a))*(!is.nan(a)) == 1]
  a = a[(!is.infinite(a))*(!is.na(a))*(!is.nan(a)) == 1]
  ta = a
  tb = ta*0
  a = c(min(a)*0.9,a,max(a)*1.1)
  b = c(0,b,1.0)
  cdf_est <- splinefun(a,b)
  
  for(i in 1:length(ta)){
    tb[i] = -qchisq(abs(0.5-cdf_est(ta[i])),df=1)*0.5 #approximate the profile 
  }
  tb = tb - min(tb)
  tb = tb/(max(tb))*max(probs)*0.4
  tb[1] = 0; tb[length(tb)] = 0; 
  polygon(ta,tb,col = alphablend(col="red",0.2),border =alphablend(col="red",0.2))  

  
  points(doses,probs,...)
  arrows(x0=doses, y0=lerror, x1=doses, 
         y1=uerror, code=3, angle=90, length=0.1)
  
}


.plot.BMDdichotomous_MA <- function(A,qprob=0.05,...){
  density_col="blueviolet"
  credint_col="azure2"
  class_list <- names(A)
  
  fit_idx    <- grep("Individual_Model",class_list)
  
  #plot the model average curve
  if ("BMDdichotomous_MA_mcmc" %in% class(A)){ # mcmc run
    n_samps <- nrow(A[[fit_idx[1]]]$mcmc_result$PARM_samples); 
    data_d   <-  A[[fit_idx[1]]]$data
    max_dose <- max(data_d[,1])
    min_dose <- min(data_d[,1])
    test_doses <- seq(min_dose,max_dose,(max_dose-min_dose)/500); 
    ma_samps <- sample(fit_idx,n_samps, replace=TRUE,prob = A$posterior_probs)
    temp_f   <- matrix(0,n_samps,length(test_doses))
    print(dim(temp_f))
    temp_bmd <- rep(0,length(test_doses))
    
    probs <- (0.5+data_d[,2,drop=T])/(1.0 + data_d[,3,drop=T])
    se <- sqrt(probs*(1-probs)/data_d[,3,drop=T])
    doses = data_d[,1,drop=T]
    uerror <- apply(cbind(probs*0+1,probs+se),1,min)
    lerror <- apply(cbind(probs*0,probs-se),1,max)
    
    dose = c(doses,doses)
    Response = c(uerror,lerror)
    plot(dose,Response,type='n',main=fit$full_model...)
    
    print(n_samps)
    print(length(fit_idx))
    for (ii in 1:n_samps){
      fit <- A[[fit_idx[ma_samps[ii]]]]
      
      if (fit$model=="hill"){
        temp_f[ii,] <- .dich_hill_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
      }
      if (fit$model=="gamma"){
        temp_f[ii,] <- .dich_gamma_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
      }
      if (fit$model == "logistic"){
        temp_f[ii,] <- .dich_logist_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
      }
      if (fit$model=="log-logistic"){
        temp_f[ii,] <- .dich_llogist_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
      }
      if (fit$model=="probit"){
        temp_f[ii,] <- .dich_probit_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]
      }
      if (fit$model=="log-probit"){
        temp_f[ii,] <- .dich_lprobit_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]   
      }
      if (fit$model=="multistage"){
        temp_f[ii,] <- .dich_multistage_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii]   
      }
      if (fit$model=="qlinear"){
        temp_f[ii,] <-  .dich_qlinear_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii] 
      }
      if (fit$model=="weibull"){
        temp_f[ii,] <- .dich_weibull_f(fit$mcmc_result$PARM_samples[ii,],test_doses)
        temp_bmd[ii] <- fit$mcmc_result$BMD_samples[ii] 
      }
    }
    me <- colMeans(temp_f)
    lq <- apply(temp_f,2,quantile, probs = qprob)
    uq <- apply(temp_f,2,quantile, probs = 1-qprob)
    col1 = alphablend(credint_col,1)
    
    polygon(c(test_doses,test_doses[length(test_doses):1]),
            c(uq,lq[length(test_doses):1]),col = col1,border=col1)
    lines(test_doses,me,lwd=2)
    temp_fit <- splinefun(test_doses,me)
    bmd <- quantile(temp_bmd,c(qprob,0.5,1-qprob),na.rm = TRUE)
    if(sum(!is.nan(test_doses) + !is.infinite(test_doses)) == 0){ 
      lines( c(bmd[1],bmd[1]),c(0,temp_fit(bmd[1])))
      lines( c(bmd[2],bmd[2]),c(0,temp_fit(bmd[2])))
      lines( c(bmd[3],bmd[3]),c(0,temp_fit(bmd[3])))
    }
    
    
    temp = temp_bmd[!is.nan(temp_bmd)]
    temp = temp[!is.infinite(temp)]
    temp = temp[temp < 30*max(doses)]
    Dens =  density(temp,cut=c(quantile(temp_bmd,0.995,na.rm = TRUE)))
    Dens$y = Dens$y/max(Dens$y) * (max(Response)-min(Response))*0.4
    temp = which(Dens$x < max(test_doses*30))
    D1_y = Dens$y[temp]
    D1_x = Dens$x[temp]
    qm = min(Response)
    polygon(c(0,D1_x,max(doses)),c(qm,qm+D1_y,qm),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
    #plot the individual models proportional to their weight
    temp_f <- rep(0,length(test_doses))
    for (ii in 1:length(fit_idx)){
      fit <- A[[fit_idx[ii]]]
      if (fit$model=="hill"){
        temp_f <- .dich_hill_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="gamma"){
        temp_f <- .dich_gamma_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model == "logistic"){
        temp_f <- .dich_logist_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="log-logistic"){
        temp_f <- .dich_llogist_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="probit"){
        temp_f <- .dich_probit_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="log-probit"){
        temp_f<- .dich_lprobit_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="multistage"){
        temp_f <- .dich_multistage_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="qlinear"){
        temp_f<- .dich_qlinear_f(fit$fitted_model$parameters,test_doses)
      }
      if (fit$model=="weibull"){
        temp_f<- .dich_weibull_f(fit$fitted_model$parameters,test_doses)
      }
      col = alphablend(col='coral3',A$posterior_probs[ii])
      lines(test_doses,temp_f,col=col,lwd = 2)
    }
    
  }else{
    
  }
  
  
  points(doses,probs)
  arrows(x0=doses, y0=lerror, x1=doses, 
         y1=uerror, code=3, angle=90, length=0.1)
  
}
