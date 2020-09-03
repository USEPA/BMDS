# FUNL
cont_FUNL_f <- function(A,doses){
     b <- A[1] + A[2]*exp((doses-A[5])^2*(-A[6]))*(1/(1+exp(-(doses-A[3])/A[4])))
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
#dichotomous log-probit
cont_exp_3_f <-function(parms,d){
  g <- parms[1]
  b <- parms[2]
  e <- parms[3] 
  rval <- g*exp((b*d)^e)
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
  
  density_col="red"
  credint_col="lightblue1"
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
  
  density_col="red"
  credint_col="lightblue1"
  
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

