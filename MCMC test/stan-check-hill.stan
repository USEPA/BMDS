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
  vector[N] y;
  vector[N] d; 
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real a;
  real b; 
  real k; 
  real n; 
  real sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a ~ lognormal(3.63653,0.1); 
  b ~ normal(0,37.96); 
  k ~ lognormal(4.07754,0.33); 
  n ~ lognormal(0,0.33);
  sigma ~ normal(0,2); 
  for (i in 1:N)
      y[i] ~ normal(a + b*pow(d[i],n)/(k^n+pow(d[i],n)), exp(0.5*sigma));
  
}

