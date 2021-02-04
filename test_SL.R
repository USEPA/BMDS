#Test Matt's code

# Continous case
doses<-c(rep(0,5),rep(25,5),rep(50,5),rep(75,5),rep(100,5))

y<-c(11.041634,9.842275, 9.184704,10.638314, 10.861213, 9.252667, 8.458513,
     10.005290,  9.835388, 8.205600, 7.547495, 8.163219, 7.510636,  8.442046,
     8.345230, 8.100379, 8.194811, 7.740399, 7.973255, 7.779894, 7.728524, 8.063397,
     7.950163, 7.655579, 8.470017)
  
#   
# model_list = data.frame(model_list = c(rep("hill",3),rep("exp-3",3),rep("exp-5",3),rep("power",2)),
#                         distribution_list =  c(rep(c("normal","normal-ncv","lognormal"),3),"normal",
#                                                "normal-ncv"))

AA <- ma_continuous_fit(as.matrix(doses),as.matrix(y), model_list = c(rep("hill",1),rep("exp-3",1),rep("exp-5",1),rep("power",1)),
                        distribution_list =  c(rep(c("normal"),3),"normal"),fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1)



cleveland_plot(AA)

plot(AA)
