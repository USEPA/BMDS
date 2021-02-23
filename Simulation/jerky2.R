library(ToxicR)
model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3"),rep("exp-5"),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),
                                                rep(c("normal","lognormal"),2),"normal",
                                                "normal-ncv"))
model_list2 = data.frame(model_list = c(rep("hill",1),rep("exp-3",1),rep("exp-5",1),rep("power",1)),
                         distribution_list =  c(rep(c("normal"),4)))

for (ii in 1:10000){
  print(ii)
  AA_l <- ma_continuous_fit(as.matrix(doses),as.matrix(y),model_list=model_list,
                            fit_type = "laplace",BMD_TYPE = 'hybrid',BMR = 0.05,point_p = 0.025)
  print(AA_l$posterior_probs)
}
q = single_continuous_fit(as.matrix(doses),as.matrix(y),model_type = "exp-5",distribution = "normal-ncv",BMD_TYPE = "hybrid",
                          BMR = 0.05,point_p = 0.025,fit_type = 'laplace',sstat = F)
q$bmd