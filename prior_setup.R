
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


a = list(hill.n,exp5.n.ncv,exp5.n)

bob = create_continuous_prior(prior,"exp-3","normal-ncv")