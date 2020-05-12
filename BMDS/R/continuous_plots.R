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

#dichotomous weibull
cont_power_f <-function(parms,d){
  g <- 1/(1+exp(-parms[1])); 
  a <- parms[2];
  b <- parms[3]; 
  rval <- g + (1-g)*(1-exp(-b*d^a))
  return (rval)
}



plot.BMDcont_fit_MCMC<-function(fit,qprob=0.05,...){
  
  density_col="red"
  credint_col="lightblue1"
  BMD_DENSITY = T
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
  
  if (ncol(fit$data) ==4 ){ #sufficient statistics
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
  if (fit$model=="hill"){
    Q <- apply(fit$PARMS,1,cont_hill_f, d=test_doses)
  }
  if (fit$model=="exp-3"){
    Q <- apply(fit$PARMS,1,cont_exp_3_f, d=test_doses)
  }
  if (fit$model=="exp-5"){
    Q <- apply(fit$PARMS,1,cont_exp_5_f, d=test_doses)
  }
  
  if (fit$model=="power"){
    Q <- apply(fit$PARMS,1,cont_power_f, d=test_doses)
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
    Dens =  density(fit$BMD,cut=c(max(test_doses)),adjust =1.5)
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
  

plot.BMDcont_fit_mle<-function(fit,qprob=0.05,...){
  
  density_col="red"
  credint_col="lightblue1"
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
   
  if (ncol(fit$data) ==4 ){ #sufficient statistics
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
  
  cdf_rows = nrow(fit$cdf)
  a = fit$cdf[,1]
  b = fit$cdf[,2]
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
  tb = tb/(max(tb))*(max(Response)-min(Response))*0.6
  tb[1] = 0; tb[length(tb)] = 0; 
  tb = tb + min(Response)
  polygon(ta,tb,col = alphablend(col="red",0.2),border =alphablend(col="red",0.2))  

  
  if (ncol(fit$data) ==4){
       points(doses,mean,...)
       arrows(x0=doses, y0=lerror, x1=doses, 
              y1=uerror, code=3, angle=90, length=0.1)
  }else{
       points(doses,Response,...)
  }
  
}




plot.BMDcont_fit_laplace<-function(fit,qprob=0.05,...){
  
  density_col="red"
  credint_col="lightblue1"
  
  
  if (qprob < 0 || qprob > 0.5){
    stop( "Quantile probability must be between 0 and 0.5")
  }
  
  
  
  if (ncol(fit$data) ==4 ){ #sufficient statistics
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
  
  cdf_rows = nrow(fit$cdf)
  a = fit$cdf[,1]
  b = fit$cdf[,2]
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
  tb = tb/(max(tb))*(max(Response)-min(Response))*0.6
  tb[1] = 0; tb[length(tb)] = 0; 
  tb = tb + min(Response)
  polygon(ta,tb,col = alphablend(col="red",0.2),border =alphablend(col="red",0.2))  
  
  
  if (ncol(fit$data) ==4){
       points(doses,mean,...)
       arrows(x0=doses, y0=lerror, x1=doses, 
              y1=uerror, code=3, angle=90, length=0.1)
  }else{
       points(doses,Response,...)
  }
  
}
