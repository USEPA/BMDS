#################################################
# Prior File
#
#
#
#################################################

normprior<-function(mean = 0, sd = 1, lb = -100,ub=100){
      if (ub < lb){
        stop("Upper Bound must be greater than lower bound")
      }
      retValue <- matrix(c(1,mean,sd,lb,ub),ncol = 5)
      class(retValue) <-"BMDprior"
      return(retValue)
}

lnormprior<-function(mean = 0, sd = 1, lb = -100,ub=100){
  if (lb < 0){
    lb = 0
  }
  
  if (ub < lb){
    stop("Upper Bound must be greater than lower bound")
  }
  
  retValue <- matrix(c(2,mean,sd,lb,ub),ncol = 5)
  class(retValue) <-"BMDprior"
  return(retValue)
}

print.BMDprior<-function(prior){
  if(prior[1] == 1){
    cat(sprintf("Prior: Normal(mu = %1.2f, sd = %1.3f) 1[%1.2f,%1.2f]\n",prior[2],
                prior[3],prior[4],prior[5]))
    return();
  }
  if (prior[1] == 2){
    cat(sprintf("Prior: Log-Normal(log-mu = %1.2f, log-sd = %1.3f) 1[%1.2f,%1.2f]\n",prior[2],
                prior[3],prior[4],prior[5]))
    return();
  }
  cat("Distribution not specified.")
}

create_prior_list <- function(x1,x2,...){
  cl <- match.call()
  mf <- as.list(match.call(expand.dots = TRUE))[-1]
  X <- matrix(0,nrow=length(mf),ncol=5)
  for (ii in 1:length(mf)){
       X[ii,] = eval(mf[[ii]])
  }
  X <- list(X)
  class(X) <- "BMDmodelprior"
  return(X)
}

combine_prior_lists<-function(p1,p2){
  if (as.character(class(p1)) == "BMDprior"){
    x1 <- as.matrix(p1[ ,,drop=F])
  
  }else{
    x1 <- p1[[1]]
  }
  
  if (as.character(class(p2)) == "BMDprior"){
    x2 <- as.matrix(p2[,,drop=F])
  
  }else{
    x2 <- p2[[1]]  
  }
  
  retval <- list(rbind(x1,x2))
  class(retval) <- "BMDmodelprior"
  return(retval)
}

print.BMDmodelprior <- function(priors){
  X = priors[[1]]
  cat("Model Parameter Priors\n ")
  cat("------------------------------------------------------------------------\n")
  for (ii in 1:nrow(X)){
    V = X[ii,]
    class(V) <- "BMDprior"
    print(V)
  }
}

#################################################33
# bayesian_prior_dich(model,variance)
##################################################
bayesian_prior_continuous  <- function(model,variance){
  
  dmodel = which(model==c("hill","exp-3","exp-5","power"))
  dvariance = which(variance == c("normal","normal-ncv","lognormal"))
  
  #Hill Prior NonConstant Normal Prior
  if (dmodel == 1 && dvariance == 2){
    prior <- create_prior_list(lnormprior(0,0.1,0,100),
                               normprior(1,2,-10,100),
                               lnormprior(log(0.5),1,0,18),
                               lnormprior(log(1.5),0.250099980007996,0,18),
                               lnormprior(0,1,0,18),
                               normprior(0, 1,-18,18));
    return(prior)
  }
  
  #Exponential NonConstant Normal Prior
  if (dmodel == 2 && dvariance == 2){
      prior <- create_prior_list(lnormprior(0,0.1,0,100),
                               lnormprior(0,1, 0,18),
                               normprior(0,1, -20,20),    # log(c)
                               lnormprior(0,0.250099980007996,0,18), #d 
                               lnormprior(0,0.5,0,18), 
                               normprior(0,1,-18,18));
      return(prior)
  }
  #Exp-5 Nonconstnat Normal 
  if (dmodel == 3 && dvariance == 2){
    prior <- create_prior_list(lnormprior(0,0.1,0,100),
                               normprior (0,1, -18,18),
                               normprior(0,1, -20,20),    # log(c)
                               lnormprior(0.405,0.250099980007996,0,18), #d 
                               lnormprior(0,0.5,0,18), 
                               normprior(0,1,-18,18));
    return(prior)
  }
  
  #Power NonConstant Normal Prior
  if (dmodel == 4 && dvariance == 2){
    prior <- create_prior_list(lnormprior(0,0.1,0,100), # a
                               normprior(0,1,  -1e4,1e4),     # b
                               lnormprior(log(1.5),0.5, 0,40),  #k
                               lnormprior(0,0.250099980007996,0,18),
                               normprior(0,1,-18,18))
    return(prior)
  }
  
  
  #Hill model
  if (dmodel == 1){
    prior <- create_prior_list(lnormprior(0,0.1,0,100),
                               normprior(0,2,-100,100),#normprior(1,2,-18,18),
                               lnormprior(log(0.5),1,0,18),
                               lnormprior(log(1.5),0.250099980007996,0,18),
                               normprior(0,1,-18,18)); 
    return(prior)
  }
  
  #Exponential 
  if (dmodel == 2){
    prior <- create_prior_list(lnormprior(0,0.1, 0,100), # a
                               lnormprior(0,1, 0,18),     # b
                               normprior(0,1, -20,20),    # log(c)
                               lnormprior(0,0.250099980007996,0,18), #d 
                               normprior(0,1,-18,18))
    return(prior)
  }
  #Power NonConstant Normal Prior
  if (dmodel == 4){
   prior <-create_prior_list(lnormprior(0,0.1,0,100), # a
                             normprior(0,1,  -1e4,1e4),     # b
                             lnormprior(log(1.5),0.5, 0,40),  #k
                             normprior(0,1,-18,18))
   return(prior)
  }
  
  #Exponential-5
  if (dmodel == 3){
    prior <- create_prior_list(lnormprior(0,0.1, 0,100), # a
                               lnormprior(0,1, 0,18),     # b
                               normprior(0,1, -20,20),    # log(c)
                               lnormprior(0,0.250099980007996,0,18), #d 
                               normprior(0,1,-18,18))
    return(prior)
  }
}

##############################################################
#Standard Dichtomous 
##############################################################
bayesian_prior_dich  <- function(model,degree=2){
  dmodel = which(model==c("hill","gamma","logistic", "log-logistic",
                          "log-probit"  ,"multistage"  ,"probit",
                          "qlinear","weibull"))
  if (dmodel==1){ #HILL
    prior <- create_prior_list(normprior(	-1,	2,	-40,	40),
                               normprior( 0,	3,	-40,	40),
                               normprior(-3,	3.3,	-40,	40),
                               lnormprior(0.693147,	0.5,	0,	40))
  }
  if (dmodel==2){ #GAMMA
    prior <- create_prior_list(normprior(	0,	2,	-18,	18),
                               lnormprior(	0.693147180559945,	0.424264068711929,	0.2,	20),
                               lnormprior(	0,	1,	0,	1e4))
    
  }
  if (dmodel == 3){ #LOGISTIC
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.1,	1,	0,	40))
  }
  if (dmodel == 4){ #LOG-LOGISTIC
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               normprior(0,	1,	-40,	40),
                               lnormprior(0.693147180559945,	0.5,	0,	20))
  }
  if (dmodel == 5){ #LOG-PROBIT
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               normprior(	0,	1,	-40,	40),
                               lnormprior(0.693147180559945,	0.5,	0,	20))
  }
  
  if (dmodel == 6){ #MULTISTAGE
    startP <- create_prior_list(normprior(	0,	2,	-20,	20),
                                lnormprior( 	0,	0.5,	0,	100))
    
    if (degree >= 2){#make sure it is a positive degree
      for (ii in (2:degree)){
        startP <- combine_prior_lists(startP,lnormprior(0,1,0,1e6))
      }
    }
    prior <- startP
  }
  if (dmodel == 7){ #PROBIT
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.1,	1,	0,	40))
  }
  if (dmodel == 8){ #QLINEAR
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.15,  1,	0,	18))
  }
  if (dmodel == 9){ #WEIBULL
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.424264068711929,	0.5,	0,	40),
                               lnormprior(0,	1.5,	0,	1e4))
  }  
  
  return(prior)
}
