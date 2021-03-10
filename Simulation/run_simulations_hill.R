library(ToxicR)
model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",3),rep("exp-5",3),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv","lognormal"),2),
                                                "normal", "normal-ncv"))
model_list2 = data.frame(model_list = c(rep("hill",1),rep("exp-3",1),rep("exp-5",1),rep("power",1)),
                         distribution_list =  c(rep(c("normal"),4)))

file_list = dir()
file_list = file_list[!(file_list %in% "results")]
options(warn=-1)
for (ii in 1:length(file_list)){
  load(file_list[ii])
  BMD_result_SD_ML1_mcmc = matrix(NA,1000,3)
  BMD_result_SD_ML1_lapl = matrix(NA,1000,3)
  BMD_result_SD_ML2_mcmc = matrix(NA,1000,3)
  BMD_result_SD_ML2_lapl = matrix(NA,1000,3)
  
  BMD_result_HB_ML1_mcmc = matrix(NA,1000,3)
  BMD_result_HB_ML1_lapl = matrix(NA,1000,3)
  BMD_result_HB_ML2_mcmc = matrix(NA,1000,3)
  BMD_result_HB_ML2_lapl = matrix(NA,1000,3)
  
  pprobs_ML1 = matrix(NA,1000,10)
  pprobs_ML2 = matrix(NA,1000,4)
  for (jj in 1:1000){#nrow(sim_data)){
      print(sprintf("File:%d Iter:%d.",ii,jj))
      ###############################################################################
      y = sim_data[jj,]
   
      AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                              fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1,samples = 75000)
      BB <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list2,
                              fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1,samples = 75000)
      AA_l <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                              fit_type = "laplace",BMD_TYPE = 'sd',BMR = 1)
      BB_l <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list2,
                              fit_type = "laplace",BMD_TYPE = 'sd',BMR = 1)
      BMD_result_SD_ML1_mcmc[jj,] = AA$bmd
      BMD_result_SD_ML1_lapl[jj,] = AA_l$bmd
      BMD_result_SD_ML2_mcmc[jj,] = BB$bmd
      BMD_result_SD_ML2_lapl[jj,] = BB_l$bmd
      
      pprobs_ML1[jj,] = AA$posterior_probs
      pprobs_ML2[jj,] = BB$posterior_probs
      ###############################################################################
      AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                              fit_type = "mcmc",BMD_TYPE = 'hybrid',BMR = 0.05,point_p = 0.025,samples = 75000)
      BB <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list2,
                              fit_type = "mcmc",BMD_TYPE = 'hybrid',BMR = 0.05,point_p = 0.025,samples = 75000)

      AA_l <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                              fit_type = "laplace",BMD_TYPE = 'hybrid',BMR = 0.05,point_p = 0.025)
      BB_l <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list2,
                              fit_type = "laplace",BMD_TYPE = 'hybrid',BMR = 0.05,point_p = 0.025)
      BMD_result_HB_ML1_mcmc[jj,] = AA$bmd
      BMD_result_HB_ML1_lapl[jj,] = AA_l$bmd
      BMD_result_HB_ML2_mcmc[jj,] = BB$bmd
      BMD_result_HB_ML2_lapl[jj,] = BB_l$bmd
     ################################################################################
    }
  
  save(BMD_result_HB_ML1_lapl,BMD_result_HB_ML2_lapl,
       BMD_result_HB_ML1_mcmc,BMD_result_HB_ML2_mcmc,
       BMD_result_SD_ML1_lapl,BMD_result_SD_ML2_lapl,
       BMD_result_SD_ML1_mcmc,BMD_result_SD_ML2_mcmc,
       pprobs_ML1,pprobs_ML2,file=sprintf("./results/simrun_%s",file_list[ii]))
}