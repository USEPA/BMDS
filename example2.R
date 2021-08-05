library(ToxicR)
library(ggplot2)

prior <- create_prior_list(normprior(0,1,-100,100),
                           normprior(0,10,-1e4,1e4),
                           lnormprior(0,1,0,100),
                           lnormprior(log(2),0.4215,0,18),
                           #lnormprior(0, 0.75,0,100),
                           normprior(0, 10,-100,100));

bob = create_continuous_prior(prior,"exp-5","lognormal")

v1 <- c(13.184152,12.8906975,12.359554,13.073001,12.861814,12.967434,12.88052,
        13.249991,	12.992931,	13.022338,	13.614057,	13.287018,	13.449239,	13.950747,
        13.239134,	13.82321,	15.080262,	14.85966,	14.7805395,	15.238369,	14.749196,
        14.913585,	15.181719,	15.051697,	15.065641,	15.16396,	15.484345,	16.493923,
        15.633442,	15.96033,	15.388061)

prior <- create_prior_list(lnormprior(0,1,-100 ,100),
                           normprior( 0, 1,-100,100),#normprior(1,2,-18,18),
                           lnormprior(0 ,1,0   ,100),
                           lnormprior(0,1 ,0   ,18 ),
                           normprior(0,2  ,-18 ,18)); 

Bob <- ma_continuous_fit(as.matrix(doses),as.matrix(v1)+10,fit_type ="laplace",BMR = 1, model_list = a )


library(readr)
PFOA_Liver <- read_table2("~/OneDrive - National Institutes of Health/Projects/2021/Continuous MA/data/PFOA_Liver_lin_Control_Ratio.txt", 
                          col_names = FALSE)


PFOA_Liver <- PFOA_Liver[-(1:2),]


doses	<- c(0,	0,	0,	0,	0.156,	0.156,	0.156,	0.3125,	0.3125,	0.3125,
           0.625,	0.625,	0.625,	1.25,	1.25,	1.25,	2.5,	2.5,	2.5,	5,5,
           5,	5,	10,	10,	10,	10,	20,	20,	20,	20) 



B <- single_continuous_fit(as.matrix(doses),as.matrix(v1)+10,model_type = "hill", fit_type = "mle")
#BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "hill", distribution="normal",fit_type = "mle",sstat = F,)

library(dplyr)
library(ToxicR)
temp <- PFOA_Liver %>% filter(X1 == "ABCG2_32656")
v1 <- as.numeric(temp[2:length(temp)])

R  <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "hill", distribution="normal",fit_type = "laplace",BMR = 1,isFast = FALSE)
R2  <- single_continuous_fit(as.matrix(doses),as.matrix(v1)/sd(v1[doses==0]),model_type = "hill", distribution="normal",fit_type = "mle",BMR = 1,isFast = FALSE)

B  <- single_continuous_fit(as.matrix(log(doses+0.001)-log(0.001)),as.matrix(v1),model_type = "FUNL", distribution="normal",fit_type = "mcmc",BMR = 1,isFast = FALSE,samples = 200000)


# I guess based on my understanding we do not need Density plot for laplace fitting case - SL 
# Aesthetically, why it doesn't show the other bmd line? 

B  <- single_continuous_fit(as.matrix(doses2),as.matrix(v1),model_type = "hill", distribution="normal-ncv",fit_type = "laplace",BMR = 2,sstat = F,isFast = FALSE)

plot(B) 

# Updated 06/02/21
.plot.BMDcont_fit_maximized(B)
B$bmd
# Multiple parameters, and they have power different cases;
# It can be repeated 

# Model average fitting 

model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))



# 05/28 SL Try to focus on the issue here
# 1. CI Band seems odd
# 2. BMD Density plot looks odd 
# 3. Color Theme update for giving user more option

# I think it is useful that which model is dominant in terms of the line- color
# Range max-min part should be updated... 
# BMDU is out of bound


# How about the log scale?
Bob <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "mcmc",BMR = 2,model_list = model_list )

# X - axis is bit odd for log10() case 
# For model average why it shows only one model?

plot(BB)

# Question not sure why BMD is superimposed? 
.plot.BMDcontinuous_MA(BB, qprob=0.05)

# X - axis is bit odd for log10() case 
.plot.BMDcontinuous_MA(BB, qprob=0.05)+scale_x_log10()



C  <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal-ncv",fit_type = "mcmc",BMR = 3)
BB <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "laplace")
AA <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),model_list=model_list,
                        fit_type = "laplace",BMR = 1,samples = 75000)
model_list  = data.frame(model_list = c(rep("hill",2),rep("exp-3",2),rep("exp-5",2),rep("power",2)),
                         distribution_list =  c(c("normal","normal-ncv"),rep(c("normal","normal-ncv"),2),
                                                "normal", "normal-ncv"))
temp <- PFOA_Liver %>% filter(X1 == "CPT1B_8373")
v1 <- as.numeric(temp[2:length(temp)])
system.time({B  <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "hill", distribution="normal",fit_type = "mle",BMR = 3,sstat = F,isFast=FALSE)})
plot(B)
g <- matrix(c( 0.272012, 0.453829, -0.125416, 0.357458,0.79383,-0.563791),nrow = 1)
sqrt(g%*%B$covariance%*%t(g))

C  <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "hill", distribution="normal",fit_type = "mle",BMR = 3,sstat = F)


BB <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "laplace",BMR=1.5)

#AA <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),fit_type = "mcmc")
#plot(BB)

temp <- PFOA_Liver %>% filter(X1 == "CYP3A1_32809")
v1 <- as.numeric(temp[2:length(temp)])
system.time({BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "hill", distribution="normal-ncv",fit_type = "mle",BMR=3,isFast = T)})
plot(BB)

g <- matrix(c(0,-0.747871,0.560669,0.703165, 0.832507), nrow = 1)
g%*%BB$covariance%*%t(g)
d <- seq(0,20,0.1)
y <- cont_hill_f(BB$parameters,d)
plot(log10(doses+0.01),v1,pch=16)
lines(log10(d+0.01),y)

v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal",fit_type = "mle",sstat = F,)
y <- cont_exp_5_f(BB$parameters,d)
plot(log10(doses+0.01),v1,pch=16)
lines(log10(d+0.01),y)

plot(BB)

v1 <- as.numeric(temp[2:length(temp)])
system.time({ single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "power", distribution="normal",fit_type = "mcmc",sstat = F,)})
plot(BB)


temp <- PFOA_Liver %>% filter(X1 == "ACADL_32776")
v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal",fit_type = "mle")

exp = cont_exp_5_f(BB$parameters,doses)
sum(dnorm(v1,exp,exp(BB$parameters[5]*0.5),log=T))

1-pchisq(-2*(BB$Deviance[1,1]-BB$Deviance[5,1]),BB$Deviance[1,2]-5)

plot(BB)

d <- seq(0,20,0.01)
y <- cont_hill_f(BB$parameters,d)
plot(log10(doses+0.01),v1,pch=16)
lines(log10(d+0.01),y)

plot(BB)



temp <- PFOA_Liver %>% filter(X1 == "LOC100912409_9106")
v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal",fit_type = "mcmc")
d <- seq(0,20,0.01)
y <- cont_hill_f(BB$parameters,d)
plot(log10(doses+0.01),v1,pch=16)
lines(log10(d+0.01),y)
plot(BB)




temp <- PFOA_Liver %>% filter(X1 == "AIF1_32886")
v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal",fit_type = "mcmc")
plot(BB)



temp <- PFOA_Liver %>% filter(X1 == "AIF1_32886")
v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal",fit_type = "mcmc")
plot(BB)

BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "power", distribution="normal",fit_type = "mle",sstat = F,)
plot(BB)

G = single_continuous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")


temp <- PFOA_Liver %>% filter(X1 == "SQLE_9935")
v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "hill", distribution="normal",fit_type = "mle",sstat = F)

y <- cont_hill_f(BB$parameters,d)
plot(log10(doses+0.01),v1,pch=16)
lines(log10(d+0.01),y)
plot(BB)

v1 <- as.numeric(temp[2:length(temp)])
BB <- single_continuous_fit(as.matrix(doses),as.matrix(v1),model_type = "exp-5", distribution="normal",fit_type = "mle",sstat = F,)

plot(BB)

mData <- matrix(c(0, 1,10,
                  0.3, 4,10,
                  1, 4,10,
                  4, 7,10),nrow=4,ncol=3,byrow=T)

library(ToxicR)
mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)

prior <- create_prior_list(normprior(	0,	10,	-18,	18),
                           lnormprior(	0.693147180559945,	0.01,	0.2,	20),
                           lnormprior(	0,	2,	0,	1e4))
#system.time({A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3], fit_type = "mcmc")})



#G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")
system.time({C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "laplace")})
system.time({C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit", fit_type = "laplace")})
system.time({C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit", fit_type = "laplace")})
system.time({C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull", fit_type = "laplace")})

system.time({B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mle")})

G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mcmc")
G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill", fit_type = "mle")


G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma", fit_type = "laplace")
G = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma", fit_type = "mle")


A = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "mcmc")

D = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "laplace")
D = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "mle")

E = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "laplace")
E = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "mle")

H = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "laplace")
H = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "mle")

I = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "laplace")
I = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "mle")

J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "qlinear",fit_type = "laplace")
J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "qlinear",fit_type = "mle")

J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit",fit_type = "laplace")
J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "probit",fit_type = "mle")

J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "laplace")
J = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "mle")



B = ma_dichotomous_fit(mData[,1],mData[,2],mData[,3],fit_type = "laplace")
plot(G)

plot(K[1,],ylim=c(0,1))
for (ii in 1:2000){
  lines(K[ii,])
}
plot(B$Individual_Model_5)


mData <- matrix(c(0, 0,5,
                  0.1, 2,5,
                  0.3, 3,5,
                  1, 4,5,
                  10,5,5),nrow=5,ncol=3,byrow=T)
                  
mData <- matrix(c(0, 0,50,
                  1, 1,50,
                  3, 2,50,
                  7, 5,50,
                  10, 20,50),nrow=5,ncol=3,byrow=T)

mData <- matrix(c(0, 0,10,
                  2, 0,10,
                  4, 0,10,
                  10, 10,10),nrow=4,ncol=3,byrow=T)

mData <- matrix(c(0, 0,50,
                  62.5, 9,50,
                  125, 8,50,
                  250, 14,50),nrow=4,ncol=3,byrow=T)


A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "gamma",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "multistage",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "logistic",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-probit",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "hill",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "weibull",fit_type = "mle")

A = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "laplace")
B = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "mcmc")
C = single_dichotomous_fit(mData[,1],mData[,2],mData[,3],model_type = "log-logistic",fit_type = "mle")
plot(B,pch=16,ylab="Probability: adenomas + carcinomas",xlab="1-Bromopropane (ppm)")
plot(C,pch=16,ylab="Probability: adenomas + carcinomas",xlab="1-Bromopropane (ppm)")

cdf_rows = nrow(C$cdf)

cdf_rows = nrow(fit$cdf)s
a = fit$cdf[,1]
b = fit$cdf[,2]
b = b[(!is.infinite(a))*(!is.na(a))*(!is.nan(a)) == 1]
a = a[(!is.infinite(a))*(!is.na(a))*(!is.nan(a)) == 1]
ta = a
tb = ta*0
a = c(min(a)*0.9,a,max(a)*1.1)
b = c(0,b,1.0)
cdf_est <- splinefun(a,b)
h = 1e-6
for(i in 1:length(ta)){
  tb[i] = (cdf_est(ta[i]+h)-cdf_est(ta[i]))/h
}
  
tb = tb/(max(tb))*0.4*0.4
dens
polygon(ta,tb,col = alphablend(col="red",0.2),border =alphablend(col="red",0.2))


lines(doses,lq)
