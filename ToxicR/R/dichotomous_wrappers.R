##########################################################
#
#
#
#
##########################################################
single_dichotomous_fit <- function(D,Y,N,model_type, fit_type = "laplace",
                                    prior="default", BMR = 0.1,
                                     alpha = 0.05, degree=2,samples = 21000,
                                     burnin = 1000){
  dmodel = which(model_type==c("hill","gamma","logistic", "log-logistic",
                                "log-probit"  ,"multistage"  ,"probit",
                                "qlinear","weibull"))
  DATA <- cbind(D,Y,N)
  o1 <- c(BMR,alpha, -9999)
  o2 <- c(1,degree)
  
  if (identical(dmodel, integer(0))){
    stop('Please specify one of the following model types: 
            "hill", "gamma", "logistic", "log-logistic"
            "log-probit", "multistage"
            "probit", "qlinear", "weibull"')
  }
  if (dmodel == 6){
    if ((o2[2] < 2) + (o2[2] > nrow(DATA)-1) > 0){
      stop('The multistage model needs to have between
               2 and nrow(DATA)-1 paramaters. If degree = 1
               use the quantal linear model.')
    }
  }
  
  fitter = which(fit_type==c("mle","laplace","mcmc"))
  if (identical(fitter, integer(0)))
  {
    stop('The fit_type variable must be either "laplace","mle", or "mcmc"\n')
  }

 
  if (fitter == 1){ #MLE fit
    bounds = bmd_default_frequentist_settings(model_type,degree)
    temp = run_single_dichotomous(dmodel,DATA,bounds,o1,o2); 
    #class(temp$bmd_dist) <- "BMD_CDF"
    te <- splinefun(temp$bmd_dist[!is.infinite(temp$bmd_dist[,1]),2],temp$bmd_dist[!is.infinite(temp$bmd_dist[,1]),1],method="hyman")
    temp$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
    temp$bounds  = bounds; 
    temp$model   = model_type; 
    temp$data    = DATA
    class(temp) <- "BMDdich_fit_maximized"
    return(temp)
  }
  
  if (fitter == 2){ #laplace fit
    prior =   bayesian_prior_dich(model_type,degree); 
    temp = run_single_dichotomous(dmodel,DATA,prior[[1]],o1,o2); 
    #class(temp$bmd_dist) <- "BMD_CDF"
    te <- splinefun(temp$bmd_dist[!is.infinite(temp$bmd_dist[,1]),2],temp$bmd_dist[!is.infinite(temp$bmd_dist[,1]),1],method="hyman")
    temp$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
    
    temp$prior = prior; 
    temp$model =  model_type; 
    temp$data = DATA
    class(temp) <- "BMDdich_fit_maximized"
    return(temp)
  }
  if (fitter ==3){
    prior =   bayesian_prior_dich(model_type,degree)
    temp = run_dichotomous_single_mcmc(dmodel,DATA[,2:3,drop=F],DATA[,1,drop=F],prior[[1]],
                                       c(BMR, alpha,samples,burnin))
    #class(temp$fitted_model$bmd_dist) <- "BMD_CDF"
    temp$options = options = c(BMR, alpha,samples,burnin) ; 
    temp$prior   = prior = list(prior = prior); 
    temp$model   = model_type; 
    temp$data    = DATA
    temp$bmd     = as.numeric(c(mean(temp$mcmc_result$BMD_samples),quantile(temp$mcmc_result$BMD_samples,c(alpha,1-alpha))))
    class(temp) <- "BMDdich_fit_MCMC"
    return(temp)
  }
 
}

.print.BMD_CDF<-function(p){
  x <- splinefun(p[!is.infinite(p[,1]),2],p[!is.infinite(p[,1]),1],method="hyman")
  cat("Approximate Quantiles for the BMD\n")
  cat("--------------------------------------------------------------\n")
  cat("1% \t 5% \t 10% \t 25% \t 50% \t 75% \t 90% \t 95% \t 99%\n")
  cat("--------------------------------------------------------------\n")

  cat(sprintf("%1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t\n",
              x(0.01),x(0.05),x(0.10),x(0.25),x(.5),x(0.75),x(0.90),x(0.95),x(0.99)))
}

print.BMDdich_fit_MCMC<-function(p){
  cat ("Benchmark Dose Estimates using MCMC. \n")
  cat (sprintf("Extra Risk: BMR:%1.2f\n",p$options[1]))
  cat (sprintf("Model Type: %s\n",p$model[1]))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  m <- mean(p$BMD)
  x <- quantile(p$BMD,c(p$options[2],1-p$options[2]))
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",m,x[1],x[2],100*(1-2*p$options[2]),"% 2-sided Confidence Interval"))
}

.print.BMDdich_fit<-function(p){
  cat ("Benchmark Dose Estimates\n")
  cat ("Approximation to the Posterior\n")
  cat (sprintf("Model Type: %s\n",p$full_model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  temp = p$bmd_dist
  temp = temp[!is.infinite(temp[,1]),]
  spfun = splinefun(temp[,2],temp[,1],method = "hyman")
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",spfun(0.5),spfun(0.05),spfun(0.95),90,"% 2-sided Confidence Interval"))
}





bmd_default_frequentist_settings <- function(model,degree=2){
  dmodel = which(model==c("hill","gamma","logistic", "log-logistic",
                          "log-probit"  ,"multistage"  ,"probit",
                          "qlinear","weibull"))
  if (dmodel==1){ #HILL
    prior <- matrix(c(0,	0,	2,	-18,	18,
                      0,	0,	2,	-18,	18,
                      0,	0,	0.5,	-18,	18,
                      0,	1,	0.250099980007996,	1.00E+00,	18),nrow=4,ncol=5,byrow=T)
  }
  if (dmodel==2){ #GAMMA
    prior <- matrix(c(0,	0,	2,	-18,	18,
                      0,	1.5,	0.424264068711929,	1,	18,
                      0,	1,	1,	0,	1000),nrow=3,ncol=5,byrow=T)
  }
  if (dmodel == 3){ #LOGISTIC
    prior <- matrix(c(0,	-2,	2,	-18,	18,
                      0,	0.1,	1,	0.00E+00,	1e4),nrow=2,ncol=5,byrow=T)
  }
  if (dmodel == 4){ #LOG-LOGISTIC
    prior <- matrix(c(0,	-1.65,	2,	-18,	18,
                      0,	-2,	1,	-18,	18,
                      0,	1.2,	0.5,	0,	1e4),nrow=3,ncol=5,byrow=T)
  }
  if (dmodel == 5){ #LOG-PROBIT
    prior <- matrix(c(0,	0.01,	2,	-18,	18,
                      0,	0,	1,	-8,	8,
                      0,	1,	0.5,	 0, 1000),nrow=3,ncol=5,byrow=T)
  }
  
  if (dmodel == 6){ #MULTISTAGE
    temp <- matrix(c(0,	-2,	2,	-18,	18,
                     0, 1,	0.25,	0,	1e4,
                     0,	0.1,	1,	0,	1.00E+06),nrow=3,ncol=5,byrow=T)
    prior <- matrix(c(0,	0.1,	1,	0,	1.00E+06),nrow=1+degree,ncol=5,byrow=T)
    prior[1:3,] <- temp; 
  }
  if (dmodel == 7){ #PROBIT
    prior <- matrix(c(0,	-2,	2,	-8,	8,
                      0,	0.1,	1,	0.00E+00,	1000),nrow=2,ncol=5,byrow=T)
  }
  if (dmodel == 8){ #QLINEAR
    prior <- matrix(c(0,	-2,	2,	-18,	18,
                      0,	1,	1,	0,	1000),nrow=2,ncol=5,byrow=T)
  }
  if (dmodel == 9){ #WEIBULL
    prior <- matrix(c(0,	0,	2,	-18,	18,
                      0,	1,	1,	1.00E+00,	50,
                      0,	1,	0.424264068711929, 1.00E-06,1000),nrow=3,ncol=5,byrow=T)
  }  
  
  return(prior)
}
# fix me - remove 

bmd_default_bayesian_prior <- function(model,degree=2){
  dmodel = which(model==c("hill","gamma","logistic", "log-logistic",
                          "log-probit"  ,"multistage"  ,"probit",
                          "qlinear","weibull"))
  if (dmodel==1){ #HILL
    prior <- matrix(c(1,	-1,	2,	-40,	40,
                      1,	 0,	3,	-40,	40,
                      1,	-3,	3.3,	-40,	40,
                      2,	0.693147,	0.5,	0,	40),nrow=4,ncol=5,byrow=T)
  }
  if (dmodel==2){ #GAMMA
    prior <- matrix(c(1,	0,	2,	-18,	18,
                      2,	0.693147180559945,	0.424264068711929,	0.2,	20,
                      2,	0,	1,	0,	1e4),nrow=3,ncol=5,byrow=T)
  }
  if (dmodel == 3){ #LOGISTIC
    prior <- matrix(c(1,	0,	2,	-20,	20,
                      2,	0.1,	1,	0,	40     ),nrow=2,ncol=5,byrow=T)
  }
  if (dmodel == 4){ #LOG-LOGISTIC
    prior <- matrix(c(1,	0,	2,	-20,	20,
                      1,	0,	1,	-40,	40,
                      2,	0.693147180559945,	0.5,	1.00E-04,	20),nrow=3,ncol=5,byrow=T)
  }
  if (dmodel == 5){ #LOG-PROBIT
    prior <- matrix(c(1,	0,	2,	-20,	20,
                      1,	0,	1,	-40,	40,
                      2,	0.693147180559945,	0.5,	1.00E-04,	40),nrow=3,ncol=5,byrow=T)
  }
  
  if (dmodel == 6){ #MULTISTAGE
    temp <- matrix(c(1,	0,	2,	-20,	20,
                     2,	0,	0.5,	1.00E-04,	100,
                     2,	0,	1,	  1.00E-04,	1.00E+06),nrow=3,ncol=5,byrow=T)
    prior <- matrix(c(2,	0,	1,	  1.00E-04,	1.00E+06),nrow=1+degree,ncol=5,byrow=T)
    prior[1:3,] <- temp; 
  }
  if (dmodel == 7){ #PROBIT
    prior <- matrix(c(1,	-2,	2,	-8,	8,
                      2,	0.1,	1,	1.00E-12,	40 ),nrow=2,ncol=5,byrow=T)
  }
  if (dmodel == 8){ #QLINEAR
    prior <- matrix(c(1,	0,  	2,-18,	18,
                      2,	0.15,  1,	0,	18),nrow=2,ncol=5,byrow=T)
  }
  if (dmodel == 9){ #WEIBULL
    prior <- matrix(c(1,	0,	2,	-20,	20,
                      2,	0.424264068711929,	0.5,	0,	40,
                      2,	0,	1.5,	0,	1e4),nrow=3,ncol=5,byrow=T)
  }  
  
  return(prior)
}