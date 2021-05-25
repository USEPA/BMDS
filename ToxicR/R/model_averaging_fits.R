#################################################
# bmd_single_continuous - Run a single BMD model
#
##################################################
ma_continuous_fit <- function(D,Y,model_list=NA, fit_type = "laplace",
                                  user_priors=NA, BMD_TYPE = "sd",
                                  BMR = 0.1, point_p = 0.01, 
                                  alpha = 0.05,samples = 21000,
                                  burnin = 1000){
  myD = Y; 
  Y = as.matrix(Y)
  D = as.matrix(D)
  
  current_models = c("hill","exp-3","exp-5","power","FUNL")
  current_dists  = c("normal","normal-ncv","lognormal")
  type_of_fit = which(fit_type == c('laplace','mcmc'))

  rt = which(BMD_TYPE==c('abs','sd','rel','hybrid'))
  if (rt == 4){
    rt = 6; 
  }
  if (identical(rt,integer(0))){
    stop("Please specify one of the following BMRF types:
    		  'abs','sd','rel','hybrid','FUNL'")
  }
  
  if (rt == 4) {rt = 6;} #internally hybrid is coded as 6	
  
  
  if (is.na(model_list[[1]][1])){
    model_list = c(rep("hill",3),rep("exp-3",3),rep("exp-5",3),rep("power",2),rep("FUNL",2))
    distribution_list = c(rep(c("normal","normal-ncv","lognormal"),3),"normal","normal-ncv","normal","normal-ncv")
    model_list = data.frame(model_list = model_list, distribution_list = distribution_list)
  }

  prior_list = ma_continuous_list(model_list$model_list,model_list$distribution_list)
  
  models <- rep(0,length(prior_list))
  dlists  <- rep(0,length(prior_list))
  priors <- list()
  permuteMat = cbind(c(1,2,3,4,5),c(6,3,5,8,10)) #c++ internal adjustment
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
    temp.fit <- lm(model_data$SSTAT[,2]~model_data$X)
  }

  #Determine if there is an increasing or decreasing trend for BMD
  is_increasing = F
  if (coefficients(temp.fit)[2] > 0){
    is_increasing = T
  }
  
  options <- c(rt,BMR,point_p,alpha, is_increasing,samples,burnin)
  
  if (fit_type == "mcmc"){
    
    temp_r <- run_continuous_ma_mcmc(priors, models, dlists,Y,D,
                                    options) 
   
    tempn <- temp_r$ma_results
    tempm <- temp_r$mcmc_runs
    #clean up the run
    idx     <- grep("Fitted_Model",names(tempn))
    temp <- list()
    jj <- 1
    for ( ii in idx){
         
         temp[[jj]] <- list()
         temp[[jj]]$mcmc_result <- tempm[[ii]]
         temp[[jj]]$fitted_model <- tempn[[ii]]
         temp[[jj]]$prior <- priors[[which(ii == idx)]]
         temp[[jj]]$data  <- cbind(D,Y)
         temp[[jj]]$model <- model_list$model_list[jj]# tolower(trimws(gsub("Model: ","",temp[[ii]]$full_model)))

         #te <- splinefun(temp[[jj]]$fitted_model$bmd_dist[!is.infinite(temp[[jj]]$fitted_model$bmd_dist[,1]),2],temp[[jj]]$fitted_model$bmd_dist[!is.infinite(temp[[jj]]$fitted_model$bmd_dist[,1]),1],method="hyman")
         data_temp = temp[[jj]]$fitted_model$bmd_dist
         data_temp = data_temp[!is.infinite(data_temp[,1]) & !is.na(data_temp[,1]),]
#         data_temp = data_temp[!is.na(data_temp[,1]),]
         temp[[jj]]$bmd     <- c(NA,NA,NA)     
         
     
         if (length(data_temp) > 0){
           ii = nrow(data_temp)
            # while(ii > 2){
            #   if (abs(data_temp[ii,2] - data_temp[ii-1,2]) < 1e-4){
            #     data_temp = data_temp[-ii,]
            #     ii = nrow(data_temp)
            #   }else{
            #     ii = ii - 1; 
            #   }
            # }
             temp[[jj]]$fitted_model$bmd_dist = data_temp
             if (length(data_temp)>10 ){
             
                  te <- splinefun(data_temp[,2,drop=F],data_temp[,1,drop=F],method="hyman")
                  temp[[jj]]$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
             }
         }
         names( temp[[jj]]$bmd ) <- c("BMD","BMDL","BMDU")
        # temp[[jj]]$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
         class(temp[[jj]]) = "BMDcont_fit_MCMC"
         jj <- jj + 1
    }
    # for (ii in idx_mcmc)

    names(temp) <- sprintf("Individual_Model_%s",1:length(priors))
   # print(tempn)
    temp$ma_bmd <- tempn$ma_bmd
   
    data_temp = temp$ma_bmd
    data_temp = data_temp[!is.infinite(data_temp[,1]) & !is.na(data_temp[,1]),]
#    data_temp = data_temp[!is.na(data_temp[,1]),]
    temp$bmd     <- c(NA,NA,NA) 
   
    if (length(data_temp) > 0){
        ii = nrow(data_temp)
        #while(ii > 2){
        #  if (abs(data_temp[ii,2] - data_temp[ii-1,2]) < 1e-4){
        #    data_temp = data_temp[-ii,]
        #    ii = nrow(data_temp)
        #  }else{
        #    ii = ii - 1; 
        #  }
        #}
        temp$ma_bmd = data_temp
        tempn$posterior_probs[is.nan(tempn$posterior_probs)] = 0
        if (length(data_temp)>10 && (abs(sum(tempn$posterior_probs) -1) <= 1e-8)){
         
         
             te <- splinefun(data_temp[,2,drop=F],data_temp[,1,drop=F],method="hyman")
             temp$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
        }else{
          #error with the posterior probabilities
          temp$bmd     <- c(Inf,0,Inf)
        }
    }
   
    names(temp$bmd) <- c("BMD","BMDL","BMDU")
    temp$posterior_probs = tempn$posterior_probs;
    class(temp) <- c("BMDcontinuous_MA","BMDcontinuous_MA_mcmc")  
    return(temp)
  }else{
  
    temp   <-  run_continuous_ma_laplace(priors, models, dlists,Y,D,
                                          options)
    
    #clean up the run
    t_names <- names(temp)
    
    idx     <- grep("Fitted_Model",t_names)
    jj <- 1
    for ( ii in idx){
         temp[[ii]]$prior <- priors[[which(ii == idx)]]
         temp[[ii]]$data  <- cbind(D,Y)
         temp[[ii]]$model <- model_list$model_list[jj]# tolower(trimws(gsub("Model: ","",temp[[ii]]$full_model)))
     
         data_temp = temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[,1]),]
         if (length(data_temp)>0){
           data_temp = data_temp[!is.na(data_temp[,1]),]
           if (nrow(data_temp)>6){
                te <- splinefun(sort(data_temp[,2,drop=F]),sort(data_temp[,1,drop=F]),method="hyman")
                temp[[ii]]$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
           }else{
                temp[[ii]]$bmd     <- c(NA,NA,NA)              
           }
         }
         names( temp[[jj]]$bmd ) <- c("BMD","BMDL","BMDU")
         names(temp)[ii] <- sprintf("Individual_Model_%s",ii)
         class(temp[[ii]]) <- "BMDcont_fit_maximized"
         jj <- jj + 1
    }
    
    temp_me <- temp$ma_bmd 
    temp$bmd <- c(NA,NA,NA)
    if (!is.null(dim(temp_me))){
        temp_me = temp_me[!is.infinite(temp_me[,1]),]
        temp_me = temp_me[!is.na(temp_me[,1]),]
        temp_me = temp_me[!is.nan(temp_me[,1]),]
        temp$posterior_probs[is.nan(temp$posterior_probs)] = 0
      if (( nrow(temp_me) > 10) && abs(sum(temp$posterior_probs) -1) <= 1e-8) 
      {
        te <- splinefun(sort(temp_me[,2,drop=F]),sort(temp_me[,1,drop=F]),method="hyman")
        temp$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
      }else{
        temp$bmd     <- c(Inf, 0, Inf)
      }
    }
    names(temp$bmd) <- c("BMD","BMDL","BMDU")
    temp$posterior_probs = temp$posterior_probs;
    class(temp) <- c("BMDcontinuous_MA","BMDcontinuous_MA_laplace") 
    return (temp)
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
      priors[[ii]] = bayesian_prior_dich(model_list[ii])
    
      model_i[ii]  = .dichotomous_model_type(model_list[ii])
    }
  
  }
  ## prepare the prior list for C++ internal 
  ## by stripping it of all 'extra' class info
  temp_priors <-list(); 
  for (ii in 1:length(priors)){
       temp_priors[ii] <- priors[[ii]]
  }
  
  data <- as.matrix(cbind(D,Y,N))
  if ( fit_type == "laplace"){
    #Laplace Run
    temp <- run_ma_dichotomous(data, temp_priors, model_i,
                               model_p, FALSE, o1, o2)
    #clean up the run
    t_names <- names(temp)
    
    idx     <- grep("Fitted_Model",t_names)
    for ( ii in idx){
         temp[[ii]]$prior <- priors[[which(ii == idx)]]
         temp[[ii]]$data  <- data
         temp[[ii]]$model <- tolower(trimws(gsub("Model: ","",temp[[ii]]$full_model)))
         if (temp[[ii]]$model =="quantal-linear" ){
              temp[[ii]]$model ="qlinear"
         }
     
         te <- splinefun(temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[,1]),2],temp[[ii]]$bmd_dist[!is.infinite(temp[[ii]]$bmd_dist[,1]),1],method="hyman")
         temp[[ii]]$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
         names(temp[ii])[1] <- sprintf("Individual_Model_%s",ii)
    }
    
    class(temp) <- c("BMDdichotomous_MA","BMDdichotomous_MA_maximized")  
  }else{
    #MCMC run
    temp_r <- run_ma_dichotomous(data, temp_priors, model_i,
                               model_p, TRUE, o1, o2)
    tempn <- temp_r$ma_results
    tempm <- temp_r$mcmc_runs
    #clean up the run
    idx     <- grep("Fitted_Model",names(tempn))
    temp <- list()
    jj <- 1
    for ( ii in idx){
         
         temp[[jj]] <- list()
         temp[[jj]]$mcmc_result <- tempm[[ii]]
         temp[[jj]]$fitted_model <- tempn[[ii]]
         temp[[jj]]$prior <- priors[[which(ii == idx)]]
         temp[[jj]]$data  <- data
         temp[[jj]]$model <- tolower(trimws(gsub("Model: ","",tempn[[ii]]$full_model)))
         temp[[jj]]$options <- c(o1,o2)
         if (temp[[jj]]$model =="quantal-linear" ){
                temp[[jj]]$model ="qlinear"
         }
         te <- splinefun(temp[[jj]]$fitted_model$bmd_dist[!is.infinite(temp[[jj]]$fitted_model$bmd_dist[,1]),2],temp[[jj]]$fitted_model$bmd_dist[!is.infinite(temp[[jj]]$fitted_model$bmd_dist[,1]),1],method="hyman")
         temp[[jj]]$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
         class(temp[[jj]]) = "BMDdich_fit_MCMC"
         jj <- jj + 1
    }
   # for (ii in idx_mcmc)
    names(temp) <- sprintf("Individual_Model_%s",1:length(temp_priors))
    temp$ma_bmd <- tempn$BMD_CDF
    te <- splinefun(temp$ma_bmd[!is.infinite(temp$ma_bmd[,1]),2],temp$ma_bmd[!is.infinite(temp$ma_bmd[,1]),1],method="hyman")
    temp$bmd   <- c(te(0.5),te(alpha),te(1-alpha))
    temp$posterior_probs = tempn$posterior_probs;
    temp$post_prob
    class(temp) <- c("BMDdichotomous_MA","BMDdichotomous_MA_mcmc")  
  }
  
  return(temp)
}
  
  
 
