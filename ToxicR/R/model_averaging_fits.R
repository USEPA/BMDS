#################################################
# bmd_single_continous - Run a single BMD model
#
##################################################
ma_continuous_fit <- function(D,Y,model_list=NA, fit_type = "laplace",
                                  user_priors=NA, BMD_TYPE = "sd",
                                  BMR = 0.1, point_p = 0.01, distribution_list = NA,
                                  alpha = 0.05,samples = 21000,
                                  burnin = 1000){
  myD = Y; 
  Y = as.matrix(Y)
  D = as.matrix(D)
  
  current_models = c("hill","exp-3","exp-5","power")
  current_dists  = c("normal","normal-ncv","lognormal")
  type_of_fit = which(fit_type == c('laplace','mcmc'))

  rt = which(BMD_TYPE==c('abs','sd','rel','hybrid'))
  if (rt == 4){
    rt = 6; 
  }
  if (identical(rt,integer(0))){
    stop("Please specify one of the following BMRF types:
    		  'abs','sd','rel','hybrid'")
  }
  
  if (rt == 4) {rt = 6;} #internally hybrid is coded as 6	
  
  
  if (is.na(model_list)){
    model_list = c(rep("hill",3),rep("exp-3",3),rep("exp-5",3),rep("power",2))
    distribution_list = c(rep(c("normal","normal-ncv","lognormal"),3),"normal","normal-ncv")
  }
  prior_list = ma_continuous_list(model_list,distribution_list)
  
  models <- rep(0,length(prior_list))
  dlists  <- rep(0,length(prior_list))
  priors <- list()
  permuteMat = cbind(c(1,2,3,4),c(6,3,5,8)) #c++ internal adjustment
  for(ii in 1:length(prior_list)){
      models[ii]   <- permuteMat[which(prior_list[[ii]]$model == current_models),2] #readjust for c++ internal
      priors[[ii]] <- prior_list[[ii]]$prior[[1]]
      dlists[ii]   <- which(prior_list[[ii]]$dist == current_dists)
  }

  
  ###################  
  DATA <- cbind(D,Y);
  if (ncol(DATA)==4){
    colnames(DATA) =  c("Dose","Resp","N","StDev")
    sstat = T
  }else if (ncol(DATA) == 2){
    colnames(DATA) =  c("Dose","Resp")
    sstat = F
  }else{
    stop("The data do not appear to be in the correct format.")
  }

  model_data = list(); 
  model_data$X = D; 
  model_data$SSTAT = DATA; 

  if (sstat == T){
    temp.fit <- lm(model_data$SSTAT[,2] ~ model_data$X,
                   weights=(1/model_data$SSTAT[,4]^2)*model_data$SSTAT[,3])
  }else{
    temp.fit <- lm(model_data$SSTAT[,1]~model_data$X)
  }

  #Determine if there is an increasing or decreasing trend for BMD
  is_increasing = F
  if (coefficients(temp.fit)[2] > 0){
    is_increasing = T
  }
  
  options <- c(rt,BMR,point_p,alpha, is_increasing,samples,burnin)
  
  if (fit_type == "mcmc"){
    
    rvals <- run_continuous_ma_mcmc(priors, models, dlists,Y,D,
                                    options) 
    rvals$ma_results$options <- options
    rvals$ma_results$data    <- DATA
    return(rvals)
  }else{
  
    rvals   <-  run_continuous_ma_laplace(priors, models, dlists,Y,D,
                                          options)
    for (ii in 1:length(prior_list)){
      rvals[[ii]]$prior = prior_list[[ii]]
    }
    class(rvals) <- "BMDcont_mafit_laplace"
  #    rvals$prior <- PR
  
  # rvals$model   <- model_type
    rvals$options <- options
    rvals$data    <- DATA
    
    return (rvals)
  }
}

.dichotomous_model_type <- function(model_name){
   # based upon the following enum in the c++ file bmdStruct.h
   # enum dich_model {d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
   #d_logprobit=5, d_multistage=6,d_probit=7,
   #d_qlinear=8,d_weibull=9};
  result = which(model_name ==c("hill","gamma","logistic","log-logistic","log-probit","multistage","probit",
                                 "qlinear","weibull"))
  if ((identical(result,integer(0)))){
    stop("The model requested to be fit is not defined.")
  }
 
  return(result)
}

#################################################
# bmd_single_continous - Run a single BMD model
#
##################################################
ma_dichotomous_fit <- function(D,Y,N,model_list=integer(0), fit_type = "laplace",
                              user_priors=integer(0), BMD_TYPE = "extra",
                              BMR = 0.1, point_p = 0.01, distribution_list = NA,
                              alpha = 0.05,samples = 21000,
                              burnin = 1000){
  
  model_p <- rep(1,9)/9; # background prior is even
  o1 <- c(BMR, alpha)
  
  if (BMD_TYPE == "extra"){
    BTYPE = 1
  }else{
    BTYPE = 2 # Added risk
  }
  
  o2 <- c(BTYPE,2,samples, burnin)
  
  priors <- list();

  if (length(model_list) < 1){

    model_list =  c("hill","gamma","logistic","log-logistic","log-probit","multistage","probit",
                    "qlinear","weibull")
    model_i = rep(0,length(model_list))
    for (ii in 1:length(model_list)){
      priors[[ii]] = bmd_default_bayesian_prior(model_list[ii])
 
      model_i[ii]  = .dichotomous_model_type(model_list[ii])
    }
  
  }
  
  data <- as.matrix(cbind(D,Y,N))
  if ( fit_type == "laplace"){
    #Laplace Run
    temp <- run_ma_dichotomous(data, priors, model_i,
                               model_p, FALSE, o1, o2)
  }else{
    #MCMC run
    temp <- run_ma_dichotomous(data, priors, model_i,
                               model_p, TRUE, o1, o2)
  }
  
  return(temp)
}
  
  
 