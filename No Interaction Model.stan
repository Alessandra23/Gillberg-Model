data {
  //Define all the data
  // number of observations
  int<lower=0> N;
  int<lower=0> Ng;
  int<lower=0> Ne;
  int genotype[Ng*Ne];
  int environment[Ng*Ne];
  // response
  vector[N] y;
  // values of Kg
  matrix[Ng,Ng] Kg;
  // priors on variance of genotype
  real sigma_g;
  // priors on variance of environment
  real sigma_e;
  // priors on sigma
  real alpha;
  real beta;
}

parameters {
  real gen[Ng];
  real env[Ne];
  vector[Ng] a_g0;
  vector<lower=0>[Ne] sigma;
}

model {
  //likelihood
  for(i in 1:N){
    y[i] ~ normal(gen[genotype[i]] + env[environment[i]], sigma);
  }
  // prioris
  sigma ~ gamma(alpha,beta);
  gen ~ normal(Kg*a_g0,sigma_g);
  a_g0 ~ normal(0, sigma_g);
  env ~ normal(0, sigma_e);
}
