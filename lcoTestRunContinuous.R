# Ensure that BMDS package is loaded.
# Uncomment the next 2 lines if pkgs are not being manually loaded
#library(Rcpp)
library(BMDS)
library(rstan)

# Input data - Continuous2.dax - Individual dose-response
M = matrix(0,nrow=26,ncol=4)
colnames(M) <- c("Dose","Resp","","")
M[,1] <- c(0,0,0,0,0,18,18,18,18,20,20,20,20,30,30,30,30,35,35,35,35,59,59,59,59,59)
M[,2] <- c(39.0,39,38.4,36.3,37.1,40.2,45.3,42.1,38.3,42.5,45.2,40.1,39.8,50.1,53.4,48.2,52.1,56.1,50.4,53.2,
           55.2,55.1,59.1,56.3,52.9,53.7)
data <- list(N=length(M[,1]),
             d = M[,1], 
             y = M[,2])

#h_fit <- stan(file="stan-check-hill.stan",data=data,
#              control = list(adapt_delta=0.9),iter=10000)

C = ma_continuous_fit(M[,1,drop=F],M[,2,drop=F],fit_type="mcmc")
Q = ma_continuous_fit(M[,1,drop=F],M[,2,drop=F])

library(BMDS)


M2           <- matrix(0,nrow=5,ncol=4)
colnames(M2) <- c("Dose","Resp","N","StDev")
M2[, 1]      <- c(0,5,20,80,200)
M2[, 2]      <- c(10.56, 10.26, 8.98, 7.56, 6.99)
M2[, 3]      <- c(25,25,25,25,25)
M2[, 4]      <- c(0.56,0.26,0.35,0.21,0.33)

C = ma_continuous_fit(D=M2[,1,drop=F],Y=M2[,2:4],BMR=1,fit_type="mcmc")
B = ma_continuous_fit(D=M2[,1,drop=F],Y=M2[,2:4],BMR=1)

c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],sstat=F,BMD_TYPE="sd",BMR=1, distribution = "normal",fit_type="laplace",model_type = "exp-5")

library(BMDS)
system.time({c = single_continuous_fit(M[,1,drop=F],M[,2,drop=F],sstat=F,BMD_TYPE="sd",BMR=1, distribution = "normal",fit_type="mcmc",model_type = "power")})
system.time({b = single_continuous_fit(M[,1,drop=F],M[,2,drop=F],sstat=F,BMD_TYPE="sd",BMR=1, distribution = "normal",fit_type="laplace",model_type = "exp-5")})
system.time({a = single_continuous_fit(M[,1,drop=F],M[,2,drop=F],sstat=F,BMD_TYPE="sd",BMR=1, distribution = "normal",fit_type="mcmc",model_type = "exp-5")})


system.time({a = single_continuous_fit(M[,1,drop=F],M[,2,drop=F],sstat=F,BMD_TYPE="sd",BMR=1, distribution = "normal",fit_type="mcmc",model_type = "hill")})

system.time({q = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1,prior = c$prior,distribution = "lognormal",fit_type="laplace",model_type = "exp-5")})


system.time({c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, distribution = "lognormal",fit_type="laplace",model_type = "exp-5")})
system.time({c1 = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, distribution = "lognormal",fit_type="mle",model_type = "exp-5")})
system.time({c2 = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, distribution = "lognormal",fit_type="mcmc",model_type = "exp-5")})


system.time({b = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="hybrid",BMR=0.1, distribution = "normal-ncv",fit_type="laplace",model_type = "hill")})
system.time({a = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="hybrid",BMR=0.1, fit_type = "mcmc",  distribution = "normal-ncv",model_type = "hill")})



# Hill Priors, Normal Bayesian, Modeled variance
priorH<-    matrix(c(2,0,1,0,18,
                     1,1,2,-18,18,
                     2,log(0.5),10,0,18,
                     2,0,10,1.00E-08,18,
                     1,0,10,-18,18,
                     2,0,1,1.00E-08,18),6,5,byrow=T)

# Exponential, Normal Bayesian, Modeled variance
priorE5<-    matrix(c(1,1,2,-1e6,1e6, # a
                      2,0,2, 0,100,     # b
                      1,0,1, -40,40,    # log(c)
                      2,0,0.250099980007996,1.00E-08,18, #d
                      1,0,0.5,-18,18,
                      2,0,0.250099980007996,1.00E-08,1),6,5,byrow=T)

# EXP 5 - Normal, 1-SD, CV
c2priorE5_1SD_opt1<-  matrix(c(0, 34.9048, 1, 0, 1e8, # a
                               0, 0.0021172, 1, 0, 1e8,     # b
                               0, 5.62604, 1, 0, 1e8,    # log(c)
                               0, 1.5301, 0, 1, 1000, #d
                               0, 2.19817 , 1, -1000, 1000) # ln-alpha
                             ,5,5,byrow=T)

# Power, Normal Bayesian, Modeled variance
priorPow<-    matrix(c(1,1,2,-1e6,1e6, # a
                       1,0,10,  0,1e4,     # b
                       2,0,3, 0,40,  #k
                       1,0,0.5,-18,18,
                       2,0,0.250099980007996,1.00E-08,1),5,5,byrow=T)

# Poly 3, constant variance (CV) Priors
# - THESE ARE JUNK VALUES. ENSURE OVERRIDE = FALSE!!!!
# - Subtract rows to decrease degree
# - Add rows to increase degree
# - Add a row for modeled variance (NCV)
# - number rows = degree + 1 + {1 if CV; 2 if NCV}
priorPoly<-    matrix(c(0,0,1,0,18,
                     0,1,2,-18,18,
                     0,log(0.5),10,0,18,
                     0,0,10,1.00E-08,18,
                     0,0,1,1.00E-08,18),5,5,byrow=T)

# Option list (index starts at 0)
# 0 - BMR type
#     - AbsoluteDev = 1,
#     - StandardDev = 2,
#     - RelativeDev  = 3,
#     - PointEstimate =4,
#     - Hybrid_Extra = 6,
#     - Hybrid_Added = 7
# 1 - BMRF
# 2 - Variance type (1= Constant variance, 2= Modeled variance)
# 3 - Lognormal distribution flag (True= lognormal, False= normal)
# 4 - Polynomial degree (only used for polynomial model)
# 5 - Direction of adversity (0=auto, 1=up, -1=down)
# 6 - Restriction
#     - not poly: 0= un-restricted, 1= restricted
#     - poly: 0= un-restricted, 1= restricted up, -1= restricted down, -9999= restricted auto
# 7 - Override built-in priors flag:
#   - TRUE: Use priors that are passed in
#   - FALSE: Ignore passed-in priors; i.e., use default model priors
# 8 - Alpha (1.0 - confidence limit)
# 9 - NOT YET IMPLEMENTED: Fixed background response (-9999= estimated)

# Restricted model options (not for poly model)
opts1 <- c(2, 1, 1, F, 0, 0, 1, F, 0.05, -9999); # Constant variance
opts2 <- c(2, 1, 2, F, 0, 0, 1, F, 0.05, -9999); # Modeled variance
# poly model options: degree = 3
opts1_poly <- c(2, 1, 1, F, 3, 0, -9999, F, 0.05, -9999); # Constant variance
# EXP-5 model: Override automatic priors
#opts1_exp5 <- c(2, 1, 1, F, 0, 0, 1, T, 0.05, -9999); # Constant variance
opts1_exp5 <- opts1

# Run Hill
# resHillOpt1D2 <- run_single_continuous(6, M, priorH, FALSE, options = opts1)
# resHillOpt1D3 <- run_single_continuous(6, M2, priorH, TRUE, options = opts1)

# Run EXP
# Continuous2.dax, Option set 1
# resExp2Opt1D2 <- run_single_continuous(2, M, priorE5, FALSE, options = opts1)
# resExp3Opt1D2 <- run_single_continuous(3, M, priorE5, FALSE, options = opts1)
# resExp4Opt1D2 <- run_single_continuous(4, M, priorE5, FALSE, options = opts1)
resExp5Opt1D2 <- run_single_continuous(5, M, c2priorE5_1SD_opt1, FALSE, options = opts1_exp5)
print(resExp5Opt1D2)
# # Continuous3.dax, Option set 2
# resExp2Opt2D2 <- run_single_continuous(2, M, priorE5, FALSE, options = opts2)
# resExp3Opt2D2 <- run_single_continuous(3, M, priorE5, FALSE, options = opts2)
# resExp4Opt2D2 <- run_single_continuous(4, M, priorE5, FALSE, options = opts2)
# resExp5Opt2D2 <- run_single_continuous(5, M, priorE5, FALSE, options = opts2)

# Polynomial model
resPoly3Opt1 <- run_single_continuous(7, M, priorPoly, FALSE, options = opts1_poly)
print(resPoly3Opt1)

