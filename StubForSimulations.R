model_list = data.frame(model_list = c(rep("hill",3),rep("exp-3",3),rep("exp-5",3),rep("power",2)),
                        distribution_list =  c(rep(c("normal","normal-ncv","lognormal"),3),"normal",
                                               "normal-ncv"))

# Toxicol Pathol. 2013 Mar; 26(1): 29–34.
#Published online 2013 Apr 22. doi: 10.1293/tox.26.29
#PMCID: PMC3620211
#PMID: 23723565
#Change Trends of Organ Weight Background Data in Sprague Dawley Rats at Different Ages
#
#
library(ToxicR)

v1 <- c(13.184152,12.8906975,12.359554,13.073001,12.861814,12.967434,12.88052,
        13.249991,	12.992931,	13.022338,	13.614057,	13.287018,	13.449239,	13.950747,
        13.239134,	13.82321,	15.080262,	14.85966,	14.7805395,	15.238369,	14.749196,
        14.913585,	15.181719,	15.051697,	15.065641,	15.16396,	15.484345,	16.493923,
        15.633442,	15.96033,	15.388061)

prior <- create_prior_list(lnormprior(0,1,-100,100),
                           normprior( 0, 1,-100,100),#normprior(1,2,-18,18),
                           lnormprior(0 ,1,0,100),
                           lnormprior(0,1,0,18),
                           normprior(0,2,-18,18)); 


doses	<- c(0,	0,	0,	0,	0.156,	0.156,	0.156,	0.3125,	0.3125,	0.3125,
           0.625,	0.625,	0.625,	1.25,	1.25,	1.25,	2.5,	2.5,	2.5,	5,5,
           5,	5,	10,	10,	10,	10,	20,	20,	20,	20) 

AA <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),model_list=model_list,
                        fit_type = "mcmc")
CC <- ma_continuous_fit(as.matrix(doses),as.matrix(v1),model_list=model_list,
                        fit_type = "laplace")
