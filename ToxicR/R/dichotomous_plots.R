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



.plot.BMDdich_fit_MCMC <-function(fit,qprob=0.05,...){
  
  density_col="red"
  credint_col="lightblue1"
  BMD_DENSITY = T
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
  

  probs <- (0.5+fit$data[,2,drop=T])/(1.0 + fit$data[,3,drop=T])
  se <- sqrt(probs*(1-probs)/fit$data[,3,drop=T])
  doses = fit$data[,1,drop=T]
  uerror <- apply(cbind(probs*0+1,probs+se),1,min)
  lerror <- apply(cbind(probs*0,probs-se),1,max)

  dose = c(doses,doses)
  Response = c(uerror,lerror)
  plot(dose,Response,type='n',main=fit$full_model...)
 
  
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
    Dens =  density(temp,cut=c(max(doses)))
    Dens$y = Dens$y/max(Dens$y) * max(probs)*0.4
    temp = which(Dens$x < max(doses))
    D1_y = Dens$y[temp]
    D1_x = Dens$x[temp]
    polygon(c(0,D1_x,max(doses)),c(0,D1_y,0),col = alphablend(col=density_col,0.2),border =alphablend(col=density_col,0.2))
  }
  points(doses,probs,...)
  arrows(x0=doses, y0=lerror, x1=doses, 
         y1=uerror, code=3, angle=90, length=0.1)
  
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
  uerror <- apply(cbind(probs*0+1,probs+se),1,min)
  lerror <- apply(cbind(probs*0,probs-se),1,max)
  
  dose = c(doses,doses)
  Response = c(uerror,lerror)
  plot(dose,Response,type='n',main=fit$full_model...)
  
  
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



