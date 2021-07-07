Prior: Normal(mu = 0.00, sd = 10.000) 1[-100.00,100.00]
Prior: Normal(mu = 0.00, sd = 10.000) 1[-10000.00,10000.00]
Prior: Log-Normal(log-mu = 0.00, log-sd = 0.500) 1[0.00,100.00]
Prior: Normal(mu = 0.50, sd = 1.000) 1[0.00,100.00]
Prior: Log-Normal(log-mu = 0.00, log-sd = 0.500) 1[0.00,100.00]
Prior: Normal(mu = 0.00, sd = 10.000) 1[-200.00,200.00]
Prior: Log-Normal(log-mu = 0.00, log-sd = 0.750) 1[0.00,18.00]
Prior: Normal(mu = 0.00, sd = 10.000) 1[-100.00,100.00]

prior <- create_prior_list(normprior(0,10,-100,100),
                           normprior(0,10,-1e4,1e4),
                           lnormprior(0,0.5,0,100),
                           lnormprior(0,0.5,0,100),
                           lnormprior(0, 0.5,0,100),
                           lnormprior(0, 10,0,100),
                           normprior(0, 2,-100,100));

prior <- create_prior_list(normprior(0,10,-100,100),
                           normprior(0,10,-1e4,1e4),
                           normprior(0,0.5,0,100),
                           lnormprior(0,0.2,0,100),
                           normprior(0, 2,-100,100));

exp5.n = create_continuous_prior(prior,"exp-5","normal")

prior <- create_prior_list(normprior(0,10,-100,100),
                           normprior(0,10,-1e4,1e4),
                           normprior(0,0.5,0,100),
                           lnormprior(0,0.2,0,100),
                           lnormprior(0,0.2,0,100),
                           normprior(0, 2,-100,100));

exp5.n.ncv = create_continuous_prior(prior,"exp-5","normal-ncv")

prior <- create_prior_list(normprior(0,10,-100,100),
                           normprior(0,10,-1e4,1e4),
                           normprior(0,0.5,0,100),
                           lnormprior(0,0.2,0,100),
                           normprior(0, 2,-100,100));

hill.n = create_continuous_prior(prior,"hill","normal")

a = list(hill.n,exp5.n.ncv,exp5.n)

prior2 <- create_prior_list(normprior(0,10,-100,100),
                            normprior(0,10,-1e4,1e4),
                            lnormprior(0,0.5,0,100),
                            lnormprior(0,0.5,0,100),
                            lnormprior(0, 0.5,0,100),
                            lnormprior(0, 10,0,100),
                           lnormprior(0,0.75,0,100),
                           normprior(0, 2,-100,100));

funl = create_continuous_prior(prior2,"FUNL","normal-ncv")


a = list(hill.n,exp5.n.ncv,exp5.n,funl)

bob = create_continuous_prior(prior,"exp-3","normal-ncv")