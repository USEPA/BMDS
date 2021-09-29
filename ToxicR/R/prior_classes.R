#################################################
# Prior File
#
#################################################
#Copyright 2020  NIEHS <matt.wheeler@nih.gov>
#   
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
#is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


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
  Y <- list()
  Y[[1]] <- X
  names(Y) <- c("priors")
  class(Y) <- "BMDmodelprior"
  return(Y)
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
  
  retval <- list(priors=rbind(x1,x2))
  
  class(retval) <- "BMDmodelprior"
  return(retval)
}

.print.BMD_Bayes_model <- function(priors){
  X = priors[[1]]
  if (!is.null(priors$model)){
    cat(priors$model," Parameter Priors\n")
  }else{
    cat("Model Parameter Priors\n ")
  }
  
  cat("------------------------------------------------------------------------\n")
  for (ii in 1:nrow(X)){
    V = X[ii,]
    if (!is.null(priors$parameters)){
      temp_text = sprintf("Prior [%s]:",priors$parameters[ii])
    }else{
      temp_text = "Prior: "
    }
    if(V[1] == 1){
      cat(sprintf("%sNormal(mu = %1.2f, sd = %1.3f) 1[%1.2f,%1.2f]\n",temp_text,V[2],
                  V[3],V[4],V[5]))
    }
    if (V[1] == 2){
      cat(sprintf("%sLog-Normal(log-mu = %1.2f, log-sd = %1.3f) 1[%1.2f,%1.2f]\n",temp_text,V[2],
                  V[3],V[4],V[5]))
    }
  
  }
}


#################################################33
# bayesian_prior_dich(model,variance)
##################################################
bayesian_prior_continuous_default <- function(model,variance,degree=2){
  
  dmodel = which(model==c("hill","exp-3","exp-5","power","FUNL","polynomial"))
  dvariance = which(variance == c("normal","normal-ncv","lognormal"))
  
  #POLYNOMIAL *BLAH*
  if (dmodel == 6){
    cat("WARNING: Polynomial models may provide unstable estimates because of
         possible non-monotone behavior.\n")
    if (dvariance == 1){
      prior <- create_prior_list(normprior(0,5,-100,100))
      
      for (ii in 1:degree){
        prior <- combine_prior_lists(prior,
                                     normprior(0,5,-100,100))
      }
      
      prior <- combine_prior_lists(prior, create_prior_list(normprior (0,1,-18,18)))
    }
    
    if (dvariance == 2){
      prior <- create_prior_list(normprior(0,5,-100,100))
      
      for (ii in 1:degree){
        prior <- combine_prior_lists(prior,
                                     normprior(0,5,-100,100))
      }
      prior <- combine_prior_lists(prior, 
                                   create_prior_list(lnormprior(0,1,0,100),
                                                     normprior (0,1,-18,18)))
      
    }
    if (dvariance == 3){
      stop("Polynomial-Log-normal models are not allowed. Please 
choose normal or normal non-constant variance. \n");
    }
    return(prior)
  }
  
  #FUNL
  if (dmodel ==5 && dvariance == 1){

         prior <- create_prior_list(normprior(0,10,-100,100),
                                    normprior(0,10,-1e4,1e4),
                                    lnormprior(0,0.5,0,100),
                                    normprior(0.5,1,0,100),
                                    lnormprior(0,0.5,0,100),
                                  normprior(0,10,-200,200),
                                  normprior (0,1,-18,18))
      return(prior)
        

  }
  
  #FUNL
  if (dmodel ==5 && dvariance == 2){

       prior <- create_prior_list(normprior(0,10,-100,100),
                                  normprior(0,10,-1e4,1e4),
                                  lnormprior(0,0.5,0,100),
                                  normprior(0.5,1,0,100),
                                  lnormprior(0,0.5,0,100),
                                  normprior(0,10,-200,200),
                                  lnormprior(0,0.5,0,18), 
                                  normprior(0,2,-18,18))
       return(prior)

  }
  #Hill Prior NonConstant Normal Prior
  if (dmodel == 1 && dvariance == 2){
    prior <- create_prior_list(normprior(1,1,-100,100),
                               normprior(0,1,-100,100),
                               lnormprior(0,2,0,100),
                               lnormprior(log(1.6),0.4214036,0,18),
                               lnormprior(0, 1,0,100),
                               normprior(-3, 1,-18,18));
    return(prior)
  }
  
  #Exponential-3 NonConstant Normal Prior
  if (dmodel == 2 && dvariance == 2){

      prior <- create_prior_list( lnormprior(0,1,-100,100),
                                  lnormprior(0,2, 0,100),
                                  lnormprior(log(1.6),0.4214036,0,18),  #d 
                                  lnormprior(0,0.5,0,18), 
                                  normprior(0,1,-30,30));
      return(prior)

  }
  #Exp-5 Nonconstnat Normal 
  if (dmodel == 3 && dvariance == 2){
    prior <- create_prior_list(lnormprior(0,1,0,100),
                               normprior(0,2, -30,30),
                               normprior(0,2, -20,20),    # log(c)
                               lnormprior(log(1.6),0.4214036,0,18), #d 
                               lnormprior(0,0.5,0,18), 
                               normprior(0,1,-30,30));
    return(prior)
  }
  
  #Power NonConstant Normal Prior
  if (dmodel == 4 && dvariance == 2){
    prior <- create_prior_list(normprior(0,1,-100,100), # a
                               normprior(0,2,  -1e4,1e4),     # b
                               lnormprior(log(1.6),0.4214036, 0,40),  #k
                               lnormprior(0,0.5,0,18),
                               normprior(0,1,-18,18))

    return(prior)
  }
  
  
  
  #Hill model
  if (dmodel == 1 && dvariance == 1){
    prior <- create_prior_list(normprior(1,1,-100,100),
                               normprior( 0, 2,-1e4,1e4),
                               lnormprior(0 ,1,0,100),
                               lnormprior(log(1.6),0.4214036,0,18),
                               normprior(0,1,-30,30)); 
    return(prior)
  }
  
  #Exponential-3
  if (dmodel == 2 && dvariance == 1){
    prior <- create_prior_list(normprior(0,1, -100,100), # a
                               lnormprior(0,2, 0,100),     # b
                               lnormprior(log(1.6),0.4214036,0,18), #d 
                               normprior(0,2,-18,18))
    return(prior)
  }
  #Power NonConstant Normal Prior
  if (dmodel == 4 && dvariance == 1){

   prior <-create_prior_list(normprior(0,1,-100,100), # a
                             normprior(0,1,  -1e2,1e2),     # b
                             lnormprior(log(1.6),0.4214036, 0,40),  #k
                             normprior(0,1,-18,18))
   return(prior)
  }
  
  #Exponential-5
  if (dmodel == 3 && dvariance == 1){
    prior <- create_prior_list(lnormprior(0,1, 0,100), # a
                               normprior(0,2, -100,100),     # b
                               normprior(0,1, -20,20),    # log(c)
                               lnormprior(log(1.6),0.4214036,0,18), #d 
                               normprior(0,1,-18,18))
    return(prior)
  }
  
  #Hill model
  if (dmodel == 1 && dvariance == 3){
    prior <- create_prior_list(normprior(1,1,-100,100),
                               normprior( 0, 2,-100,100),#normprior(1,2,-18,18),
                               lnormprior(0 ,1,0,100),
                               lnormprior(log(1.6),0.4214036,0,18),
                               normprior(0,1,-18,18)); 
    return(prior)
  }
  
  #Exponential-3 
  if (dmodel == 2 && dvariance == 3){
    prior <- create_prior_list(lnormprior(0,1, 0,100), # a
                               lnormprior(0,1, 0,100),     # b
                               lnormprior(log(1.6),0.4214036,0,18), #d 
                               normprior(0,1,-18,18))

    return(prior)
  }
  
  
  #Exponential-5
  if (dmodel == 3 && dvariance == 3){
    prior <- create_prior_list(lnormprior(0,1, 0,100), # a
                               lnormprior(0,2, 0,100),     # b
                               normprior(0,2, -20,20),    # log(c)
                               lnormprior(log(1.6),0.6,0,18), #d 
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
    prior <- create_dichotomous_prior(prior,"hill")
  }
  if (dmodel==2){ #GAMMA
    prior <- create_prior_list(normprior(	0,	2,	-18,	18),
                               lnormprior(	0.693147180559945,	0.424264068711929,	0.2,	20),
                               lnormprior(	0,	1,	0,	1e4))
    prior <- create_dichotomous_prior(prior,"gamma")
  }
  if (dmodel == 3){ #LOGISTIC
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.1,	1,	0,	40))
    prior <- create_dichotomous_prior(prior,"logistic")
  }
  if (dmodel == 4){ #LOG-LOGISTIC
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               normprior(0,	1,	-40,	40),
                               lnormprior(0.693147180559945,	0.5,	0,	20))
    prior <- create_dichotomous_prior(prior,"log-logistic")
  }
  if (dmodel == 5){ #LOG-PROBIT
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               normprior(	0,	1,	-40,	40),
                               lnormprior(0.693147180559945,	0.5,	0,	20))
  }
  
  if (dmodel == 6){ #MULTISTAGE
     startP <- create_prior_list(normprior(	0,	2,	-20,	20),
                                lnormprior( 	0,	0.5,	0,	100))
    degree = floor(degree)
    if (degree >= 2){#make sure it is a positive degree
      for (ii in (2:degree)){
        startP <- combine_prior_lists(startP,lnormprior(0,1,0,1e6))
      }
    }
    prior <- startP
    prior <- create_dichotomous_prior(prior,"multistage")
  }
  if (dmodel == 7){ #PROBIT
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.1,	1,	0,	40))
    prior <- create_dichotomous_prior(prior,"probit")
  }
  if (dmodel == 8){ #QLINEAR
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.15,  1,	0,	18))
    prior <- create_dichotomous_prior(prior,"qlinear")
  }
  if (dmodel == 9){ #WEIBULL
    prior <- create_prior_list(normprior(	0,	2,	-20,	20),
                               lnormprior(0.424264068711929,	0.5,	0,	40),
                               lnormprior(0,	1.5,	0,	1e4))
    prior <- create_dichotomous_prior(prior,"weibull")
  }  
  
  return(prior)
}

#################################################33
# bayesian_prior_dich(model,variance)
##################################################
MLE_bounds_continuous  <- function(model,variance,degree=2){
  
  dmodel = which(model==c("hill","exp-3","exp-5","power","FUNL","polynomial"))
  dvariance = which(variance == c("normal","normal-ncv","lognormal"))
  
  #POLYNOMIAL *BLAH*
  if (dmodel == 6){

    if (dvariance == 1){
      prior <- create_prior_list(normprior(0,5,-1000,1000))
      
      for (ii in 1:degree){
        prior <- combine_prior_lists(prior, normprior(0,5,-18,0))
      }
      
      prior <- combine_prior_lists(prior, create_prior_list(normprior (0,1,-18,18)))
      prior[[1]][,1] = 0
      return(prior)
    }
    
    if (dvariance == 2){
      prior <- create_prior_list(normprior(0,5,0,1000))
      
      for (ii in 1:degree){
        prior <- combine_prior_lists(prior,normprior(0,5,-10000,10000))
      }
      prior <- combine_prior_lists(prior, 
                                   create_prior_list(lnormprior(0,1,0,100),
                                                     normprior (0,1,-18,18)))
      prior[[1]][,1] = 0
      return(prior)
      
    }
    if (dvariance == 3){
      stop("Polynomial-Log-normal models are not allowed. Please 
            choose normal or normal non-constant variance.\n");
    }
    return(prior)
  }
  
  #FUNL
  if (dmodel ==5 ){
    stop("The FUNL model is not available for MLE.\n");
  }
  
  
  #Hill Prior NonConstant Normal Prior
  if (dmodel == 1 && dvariance == 2){
    prior <- create_prior_list(normprior(0,1,-100,100),
                               normprior(0,2,-100,0),
                               normprior(0,1,0,18),#lnormprior(0,1,0,18),
                               lnormprior(log(1.2),1,1,18),
                               lnormprior(-2, 1,0,100),
                               normprior(0, 2,-18,18));
  }
  
  #Exponential NonConstant Normal Prior
  if (dmodel == 2 && dvariance == 2){
    prior <- create_prior_list(normprior(0,1,-0,100),
                               lnormprior(0,0.5, 0,100),
                               normprior(0,0.5, -20,20),    # log(c)
                               lnormprior(0,0.3,1,18),  #d 
                               lnormprior(0,0.5,0,18), 
                               normprior(0,1,-18,18));

  }
  #Exp-5 Nonconstnat Normal 
  if (dmodel == 3 && dvariance == 2){
    prior <- create_prior_list(normprior(0,0.1,-100,100),
                               normprior(0,1, 0,100),
                               normprior(0, 0.5, -20,20),    # log(c)
                               lnormprior(0,0.2,1,18), #d 
                               lnormprior(0,0.5,0,18), 
                               normprior(0,1,-18,18));
  
  }
  
  #Power NonConstant Normal Prior
  if (dmodel == 4 && dvariance == 2){
    prior <- create_prior_list(normprior(0,1,-100,100), # a
                               normprior(0,1,  -1e4,1e4),     # b
                               lnormprior(0,0.2,1,18),  #k
                               lnormprior(0,0.250099980007996,0,18),
                               normprior(0,1,-18,18))
   
  }
  
  #Hill model
  if (dmodel == 1 &&  dvariance == 1){
    prior <- create_prior_list(normprior(0,1,-100,100),
                               normprior(0,2,-100,0),#normprior(1,2,-18,18),
                               lnormprior(0 ,1,0,18),
                               lnormprior(1,1.2,1,18),
                               normprior(0,1,-18,18)); 
  }
  
  #Exponential 
  if (dmodel == 2 && dvariance == 1){
    prior <- create_prior_list(normprior(0,0.1, 0,100), # a
                               lnormprior(0,1, 0,100),     # b
                               normprior(0,0.5, -20,20),    # log(c)
                               lnormprior(1,0.2,1,18), #d 
                               normprior(0,1,-18,18))
  }
  #Power  Normal Prior
  if (dmodel == 4 && dvariance == 1){
    prior <-create_prior_list(normprior(0,0.1,-100,100), # a
                              normprior(0,1,  -1e2,1e2),     # b
                              lnormprior(1,0.2,0,18),  #k
                              normprior(0,1,-18,18))

  }
  
  #Exponential-5
  if (dmodel == 3 && dvariance == 1){
    prior <- create_prior_list(normprior(0,0.1, -100,100), # a
                               lnormprior(0,1, 0,100),     # b
                               normprior(0,1, -10,10),    # log(c)
                               lnormprior(1,0.2,0,18), #d 
                               normprior(0,1,-18,18))
  }
  
  #Exponential-5 lognormal
  if (dmodel == 3 && dvariance == 3){
    prior <- create_prior_list(normprior(0,0.1, -100,100), # a
                               lnormprior(0,1, 0,100),     # b
                               normprior(0,1, -10,10),    # log(c)
                               lnormprior(1,0.2,1,18), #d 
                               normprior(0,1,-18,18))
  }
  
  #Exponential-5 lognormal
  if (dmodel == 2 && dvariance == 3){
    prior <- create_prior_list(normprior(0,0.1, -100,100), # a
                               lnormprior(0,1, 0,100),     # b
                               normprior(0,1, -20,20),    # log(c)
                               lnormprior(1,0.2,1,18), #d 
                               normprior(0,1,-18,18))
  }
  
  #Hill model lognormal
  if (dmodel == 1 &&  dvariance == 3){
    prior <- create_prior_list(normprior(0,1,-100,100),
                               normprior( 0, 1,-100,100),#normprior(1,2,-18,18),
                               lnormprior(0 ,1,0,18),
                               lnormprior(1,0.2,0,18),
                               normprior(0,1,-18,18)); 
  }
  prior <- unclass(prior)
  prior[[1]][,1] = 0
  return(prior)
}
