#################################################33
#get_prior_list <- a list of Default priors for an analysis
#
#
#
##################################################


#################################################
# bmd_single_continous - Run a single BMD model
#
##################################################
single_continuous_fit <- function(D,Y,model_type="hill", fit_type = "laplace",
                                   prior="default", BMD_TYPE = "sd", sstat = T,
                                   BMR = 0.1, point_p = 0.01, distribution = "normal-ncv",
                                   alpha = 0.05,samples = 51000,
                                   burnin = 1000){
    myD = Y; 
    type_of_fit = which(fit_type == c('laplace','mle','mcmc'))
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
    dis_type = which(distribution  == c("normal","normal-ncv","lognormal"))
    dmodel = which(model_type==c("hill","exp-3","exp-5","power"))

    if (identical(dmodel, integer(0))){
      stop('Please specify one of the following model types: \n
            "hill","exp-3","exp-5","power"')
    }
    
    if (identical(rt,integer(0))){
      stop('Please specify one of the following Benchmark Dose types \n
            "hill","exp-5","power"')
    }
    
    if(identical(dis_type,integer(0))){
      stop('Please specify the distribution as one of the following:\n
            "normal","normal-ncv","lognormal"')
    }
    
    DATA <- cbind(D,Y); 
    if (ncol(DATA)==4){
      colnames(DATA) =  c("Dose","Resp","N","StDev")
    }else if (ncol(DATA) == 2){
      colnames(DATA) =  c("Dose","Resp")
    }else{
      stop("The data do not appear to be in the correct format.")
    }
    
    permuteMat = cbind(c(1,2,3,4),c(6,3,5,8))
    fitmodel = permuteMat[dmodel,2]
    
    if (prior =="default"){
      PR = bayesian_prior_continuous(model_type,distribution)
    }else{
      PR = prior; 
    }
    
    if (type_of_fit == 2){
      PR[[1]][,1] = 0
      PR[[1]][,2] = 1
    }
    
    if (distribution == "lognormal"){
      is_log_normal = TRUE
    }else{
      is_log_normal = FALSE
    }
    
    if (distribution == "normal-ncv"){
      constVar = FALSE
      vt = 2;
    }else{
      vt = 1
      constVar = TRUE
    }
    

    if (identical(rt,integer(0))){
    	 stop("Please specify one of the following BMRF types:
    		  'abs','sd','rel','hybrid'")
    }
    
    if (rt == 4) {rt = 6;} #internally hybrid is coded as 6	
    
    
    if (identical(dmodel, integer(0))){
      stop('Please specify one of the following model types: 
            "hill", "exp-3", "exp-5", "power"')
    }
    
 
    
    #check for sufficient statistics and normalize
    #if (sstat==T){
    model_data = list(); 
    model_data$X = D; model_data$SSTAT = Y
    #}else{
    #  model_data = createSuffStat(DATA[,1],DATA[,2],is_log_normal)	
    #}
    
    #start here tomorrow <- worry about this change. 
 
    if (sstat == T){
      temp.fit <- lm(model_data$SSTAT[,1] ~ model_data$X,
  		                weights=(1/model_data$SSTAT[,3]^2)*model_data$SSTAT[,2])
    }else{
      temp.fit <- lm(model_data$SSTAT[,1]~model_data$X)
    }

    #Determine if there is an increasing or decreasing trend for BMD
    is_increasing = F
    
    if (coefficients(temp.fit)[2] > 0){
	     is_increasing = T
    }
    tmodel <- permuteMat[dmodel,2] 
    
    options <- c(rt,BMR,point_p,alpha, is_increasing,constVar,point_p,samples)
    
    ## pick a distribution type 
    if (is_log_normal){
      dist_type = 3 #log normal
    }else{
      if (constVar){
          dist_type = 1 # normal
      }else{
          dist_type = 2 # normal-ncv
      }
    }
    
    if (fit_type == "mcmc"){
      
      rvals <- run_continuous_single_mcmc(fitmodel,model_data$SSTAT,model_data$X,
                                          PR[[1]],options, is_log_normal, sstat) 
      if (model_type == "exp-3"){
        rvals$PARMS = rvals$PARMS[,-3]
      }
      class(rvals) <- "BMDcont_fit_MCMC"; 
      rvals$model <- model_type
      rvals$options <- options
      rvals$data <- DATA
      rvals$bmd <- c(mean(rvals$BMD,na.rm=TRUE),quantile(rvals$BMD,c(alpha,1-alpha),na.rm=TRUE))
      rvals$prior <- PR
      return(rvals)
    }else{
      
      
      rvals   <- run_continuous_single(fitmodel,model_data$SSTAT,model_data$X,
  						                          PR[[1]],options, dist_type)
      if (type_of_fit == 2){
        #MLE was chosen
        class(rvals) <- "BMDcont_fit_mle"
        
      }else{
        #Laplace was chosen
        class(rvals) <- "BMDcont_fit_laplace"
        rvals$prior <- PR
      }
      rvals$model   <- model_type
      rvals$options <- options
      rvals$data    <- DATA
      
      return (rvals)
    }
}

print.BMDcont_fit_MCMC<-function(p){
  
  BMDtype <- c('Absolute Deviation','Standard Deviation','Relative Deviation','Hybrid')
  
  
  cat ("Benchmark Dose Estimates using MCMC. \n")
  cat (sprintf("Continuous %s BMD: BMRF-%1.2f\n",BMDtype[p$options[1]],p$options[2]))
  cat (sprintf("Model Type: %s\n",p$model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  m <- mean(p$BMD)
  x <- quantile(p$BMD,c(p$options[4],1-p$options[4]))
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",m,x[1],x[2],100*(1-2*p$options[4]),"% 2-sided Confidence Interval"))
}

print.BMDcont_fit_laplace<-function(p){
  BMDtype <- c('Absolute Deviation','Standard Deviation','Relative Deviation','Hybrid')
  
  cat ("Benchmark Dose Estimates using Laplace \n")
  cat ("approximation to the Posterior\n")
  cat (sprintf("Continuous %s BMD: BMRF-%1.2f\n",BMDtype[p$options[1]],p$options[2]))
  cat (sprintf("Model Type: %s\n",p$model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",p$bmd[1],p$bmd[2],p$bmd[3],100*(1-2*p$options[4]),"% 2-sided Confidence Interval"))
}

print.BMDcont_fit_mle<-function(p){
  BMDtype <- c('Absolute Deviation','Standard Deviation','Relative Deviation','Hybrid')
  
  cat ("Benchmark Dose Estimates using MLE \n")
  cat (sprintf("Continuous %s BMD: BMRF-%1.2f\n",BMDtype[p$options[1]],p$options[2]))
  
  cat (sprintf("Model Type: %s\n",p$model))
  cat ("BMD  (BMDL,BMDU) \n")
  cat ("---------------------\n")
  cat (sprintf("%1.2f (%1.2f,%1.2f)\n%1.2f%s\n",p$bmd[1],p$bmd[2],p$bmd[3],100*(1-2*p$options[4]),"% 2-sided Confidence Interval"))
}


#####################################################
# bmd_ma - runs the BMD continous model averaging code
#
#
####################################################
bmd_ma_continuous <- function(riskType,DATA,PR=NA,sstat=F,BMRF = 1,
							                bkg_prob=0.01,alpha=0.05){
    
    rt = which(riskType==c('AbsDev','STDev','RelDev','Point','Hybrid'))
    if (rt == 5) {rt = 6;} #internally hybrid is coded as 6	
    
    dose_scale = max(DATA[,1])
    DATA[,1] = DATA[,1]/dose_scale;   

    if (identical(rt,integer(0))){
      	 stop("Please specify one of the following BMRF types:
        		  'AbsDev','STDev','RelDev','Point','Hybrid'")
    }
    #check for sufficient statistics and normalize
    if (sstat==T){
		      model_data =    cleanSuffStat(DATA,F)
		      model_dataLN  = cleanSuffStat(DATA,T)
    }else{
		      model_data = createSuffStat(DATA[,1],DATA[,2],F)
		      model_dataLN = createSuffStat(DATA[,1],DATA[,2],T)
    }
    #Determine if there is an increasing or decreasing trend for BMD
    #estimation
    temp.fit <- lm(model_data$SSTAT[,1] ~ model_data$X,
		weights=(1/model_data$SSTAT[,3]^2)*model_data$SSTAT[,2])
 
    is_increasing = F
    if (coefficients(temp.fit)[2] > 0){
	    is_increasing = T
    }
    
    #########################################################################   
  	if (is.na(PR)) { 
  		PR = get_prior_list()
  	} else{
  			#to do add some checking to make sure the priors work
  	}
   #Options are: 
   #rt - risk type, BMRF - Benchmark Response function, 
   #bkg_prob - background probability , only used in hybrid
   #alpha - CL alpha for BMD
   #is_increasing determine if the data are increasing
   options <- c(rt,BMRF,bkg_prob,alpha,is_increasing,1);
   rvals   <- run_continuous_ma(PR, model_data$SSTAT, model_dataLN$SSTAT,model_data$X,
                								matrix(1/15,nrow=15,ncol=1), options)
   rvals$estimates$`Exp3  Normal-NCV` <- rvals$estimates$`Exp3  Normal-NCV`[-3,,drop=F]; 
   rvals$covariance$`Exp3  Normal-NCV`<-rvals$covariance$`Exp3  Normal-NCV`[-3,,drop=F]
   rvals$covariance$`Exp3  Normal-NCV`<-rvals$covariance$`Exp3  Normal-NCV`[,-3,drop=F]
   rvals$estimates$`Exp3  Normal-CV` <- rvals$estimates$`Exp3  Normal-CV`[-3,,drop=F]; 
   rvals$covariance$`Exp3  Normal-CV`<-rvals$covariance$`Exp3  Normal-CV`[-3,,drop=F]
   rvals$covariance$`Exp3  Normal-CV`<-rvals$covariance$`Exp3  Normal-CV`[,-3,drop=F]
   rvals$estimates$`Exp3  Log-Normal` <- rvals$estimates$`Exp3  Log-Normal`[-3,,drop=F]; 
   rvals$covariance$`Exp3  Log-Normal`<-rvals$covariance$`Exp3  Log-Normal`[-3,,drop=F]
   rvals$covariance$`Exp3  Log-Normal`<-rvals$covariance$`Exp3  Log-Normal`[,-3,drop=F]
   
   SCALE = model_data$SCALE; 
   temp <- clean_parameters('hill',rvals[[1]][[1]],rvals[[2]][[1]],SCALE,dose_scale,2,F)
   rvals[[1]][[1]] <- temp[[1]]; rvals[[2]][[1]] <- temp[[2]];
   temp <- clean_parameters('hill',rvals[[1]][[2]],rvals[[2]][[2]],SCALE,dose_scale,2,F)
   rvals[[1]][[2]] <- temp[[1]]; rvals[[2]][[2]] <- temp[[2]]; 
   temp <- clean_parameters('hill',rvals[[1]][[3]],rvals[[2]][[3]],SCALE,dose_scale,2,T)
   rvals[[1]][[3]] <- temp[[1]]; rvals[[2]][[3]] <- temp[[2]]; 
   temp <- clean_parameters('power',rvals[[1]][[4]],rvals[[2]][[4]],SCALE,dose_scale,2,F)
   rvals[[1]][[4]] <- temp[[1]]; rvals[[2]][[4]] <- temp[[2]]; 
   temp <- clean_parameters('power',rvals[[1]][[5]],rvals[[2]][[5]],SCALE,dose_scale,2,F)
   rvals[[1]][[5]] <- temp[[1]]; rvals[[2]][[5]] <- temp[[2]]; 
   temp <- clean_parameters('power',rvals[[1]][[6]],rvals[[2]][[6]],SCALE,dose_scale,2,T)
   rvals[[1]][[6]] <- temp[[1]]; rvals[[2]][[6]] <- temp[[2]]; 
   temp <- clean_parameters('exp-3',rvals[[1]][[7]],rvals[[2]][[7]],SCALE,dose_scale,2,F)
   rvals[[1]][[7]] <- temp[[1]]; rvals[[2]][[7]] <- temp[[2]]; 
   temp <- clean_parameters('exp-3',rvals[[1]][[8]],rvals[[2]][[8]],SCALE,dose_scale,2,F)
   rvals[[1]][[8]] <- temp[[1]]; rvals[[2]][[8]] <- temp[[2]]; 
   temp <- clean_parameters('exp-3',rvals[[1]][[9]],rvals[[2]][[9]],SCALE,dose_scale,2,T)
   rvals[[1]][[9]] <- temp[[1]]; rvals[[2]][[9]] <- temp[[2]]; 
   temp <- clean_parameters('exp-5',rvals[[1]][[10]],rvals[[2]][[10]],SCALE,dose_scale,2,F)
   rvals[[1]][[10]] <- temp[[1]]; rvals[[2]][[10]] <- temp[[2]]; 
   temp <- clean_parameters('exp-5',rvals[[1]][[11]],rvals[[2]][[11]],SCALE,dose_scale,2,F)
   rvals[[1]][[11]] <- temp[[1]]; rvals[[2]][[11]] <- temp[[2]]; 
   temp <- clean_parameters('exp-5',rvals[[1]][[12]],rvals[[2]][[12]],SCALE,dose_scale,2,T)
   rvals[[1]][[12]] <- temp[[1]]; rvals[[2]][[12]] <- temp[[2]]; 
   temp <- clean_parameters('poly',rvals[[1]][[13]],rvals[[2]][[13]],SCALE,dose_scale,2,F)
   rvals[[1]][[13]] <- temp[[1]]; rvals[[2]][[13]] <- temp[[2]]; 
   temp <- clean_parameters('poly',rvals[[1]][[14]],rvals[[2]][[14]],SCALE,dose_scale,2,F)
   rvals[[1]][[14]] <- temp[[1]]; rvals[[2]][[14]] <- temp[[2]]; 
   temp <- clean_parameters('poly',rvals[[1]][[15]],rvals[[2]][[15]],SCALE,dose_scale,2,T)
   rvals[[1]][[15]] <- temp[[1]]; rvals[[2]][[15]] <- temp[[2]]; 

   rvals$dose_scale = dose_scale; 
   rvals$SCALE      = model_data$SCALE; 
   
   return (rvals)
}

