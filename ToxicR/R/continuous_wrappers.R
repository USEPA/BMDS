#################################################33
#get_prior_list <- a list of Default priors for an analysis
#
##################################################

#################################################
# bmd_single_continous - Run a single BMD model
#
##################################################
single_continuous_fit <- function(D,Y,model_type="hill", fit_type = "laplace",
                                   prior="default", BMD_TYPE = "sd", sstat = T,
                                   BMR = 0.1, point_p = 0.01, distribution = "normal-ncv",
                                   alpha = 0.05,samples = 21000,
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
    dmodel = which(model_type==c("hill","exp-3","exp-5","power","FUNL"))

    if (identical(dmodel, integer(0))){
      stop('Please specify one of the following model types: \n
            "hill","exp-3","exp-5","power","FUNL"')
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
    #permute the matrix to the internal C values
    # Hill = 6, Exp3 = 3, Exp5 = 5, Power = 8, 
    # FUNL = 10
    permuteMat = cbind(c(1,2,3,4,5),c(6,3,5,8,10))
    fitmodel = permuteMat[dmodel,2]
    if (fitmodel == 10 && dis_type == 3){
         stop("The FUNL model is currently not defined for Log-Normal distribution.")
    }
    
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
    
    
    #Temporary fit to determine direction. 
    model_data = list(); 
    model_data$X = D; model_data$SSTAT = Y
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
      rvals$bmd <- c(mean(rvals$mcmc_result$BMD_samples,na.rm=TRUE),quantile(rvals$mcmc_result$BMD_samples,c(alpha,1-alpha),na.rm=TRUE))
      names(rvals$bmd) <- c("BMD","BMDL","BMDU")
      rvals$prior <- PR
      class(rvals) <- "BMDcont_fit_MCMC"
      return(rvals)
    }else{
      
      
      rvals   <- run_continuous_single(fitmodel,model_data$SSTAT,model_data$X,
  						                          PR[[1]],options, dist_type)
      
      rvals$bmd_dist = rvals$bmd_dist[!is.infinite(rvals$bmd_dist[,1]),]
      te <- splinefun(rvals$bmd_dist[,2],rvals$bmd_dist[,1],method="hyman")
    
      rvals$bmd     <- c(te(0.5),te(alpha),te(1-alpha))
      names(rvals$bmd) <- c("BMD","BMDL","BMDU") 
      rvals$model   <- model_type
      rvals$options <- options
      rvals$data    <- DATA
      class(rvals)  <- "BMDcont_fit_maximized"
      
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



