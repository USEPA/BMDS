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
  #define CONTINUOUS_BMD_ABSOLUTE     1
  #define CONTINUOUS_BMD_STD_DEV      2
  #define CONTINUOUS_BMD_REL_DEV      3
  #define CONTINUOUS_BMD_POINT        4
  #define CONTINUOUS_BMD_EXTRA        5
  #define CONTINUOUS_BMD_HYBRID_EXTRA 6
  #define CONTINUOUS_BMD_HYBRID_ADDED 7
  #define CONTINUOUS_BMD_EMPTY        0.0
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
 # print(D,Y)
  DATA <- cbind(D,Y);
#  print("Now I am here.s")
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
#    if (model_type == "exp-3"){
#      rvals$PARMS = rvals$PARMS[,-3]
#    }
#    class(rvals) <- "BMDcont_fit_MCMC"; 
#    rvals$model <- model_type
#    rvals$options <- options
#    rvals$data <- DATA
#    rvals$bmd <- c(mean(rvals$BMD,na.rm=TRUE),quantile(rvals$BMD,c(alpha,1-alpha),na.rm=TRUE))
  #  rvals$prior <- PR
 #   return(rvals)
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
