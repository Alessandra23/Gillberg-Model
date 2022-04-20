data {
  //Define all the data
  // number of observations
  int<lower=0> N;
  int<lower=0> I;
  int<lower=0> J;
  int aa[N];
  int bb[N];
  // response
  real y[N];
}


parameters {
  real alpha[I];
  real beta[J];
}

model {
  //likelihood
  for(i in 1:N){
    real mu = alpha[aa[i]] + beta[bb[i]];
    y[i] ~ normal(mu, 10);
  }
}
