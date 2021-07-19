//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  real scale;
  vector[N] y;
  vector[N] x; 
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real a;
  real b;
  real c;
  real<lower=0> d; 
  real<lower=0> sigma;
}
//prior <- create_prior_list(lnormprior(0,1,-100,100),
//                           normprior( 0, 1,-100,100),#normprior(1,2,-18,18),
//                           lnormprior(0 ,1,0,100),
//                           lnormprior(0,1,0,18),
//                           normprior(0,2,-18,18)); 
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a     ~ lognormal(0,1+log(scale));
  b     ~ normal(0*scale,scale*1); 
  c     ~ lognormal(0,1); 
  c     ~ lognormal(0,1); 
  sigma ~ normal(2*log(scale),1); 
  
  for ( i in 1:N){
     y[i] ~ normal(a + b*pow(x[i],d)/(pow(c,d) + pow(x[i],d)), exp(2*sigma));
  }
}

