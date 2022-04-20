data {
  //Define all the data
  // number of observations
  int<lower=0> N;
  int<lower=0> Ng;
  int<lower=0> Ne;
  int R;
  // response
  vector[N] y;
  // index of genotype
  int genotype[Ng*Ne];
  // index of environment
  int environment[Ng*Ne];
  // values of Kg
  matrix[N,R] Kg;
  // values of Ke
  matrix[N,R] Ke;
  // priors on variance of genotype
  real sigma_g;
  // priors on variance of environment
  real sigma_e;
  // priors on kernel regression weight of genotype
  real sigma_ag;
 // priors on kernel regression weight of environment
  real sigma_ae;
  // priors on sigma
  real alpha;
  real beta;
}


parameters {
  vector[Ng] gen;
  vector[Ne] env;
  vector[Ng] a_g0;
  matrix[Ng, R] ag;
  matrix[Ne, R] ae;
  matrix[Ng, R] hg;
  matrix[Ne, R] he;
  vector<lower=0>[Ne] sigma;
}

transformed parameters{
  matrix[Ng, Ne] xi;
  xi = hg * he';
}

model {
   //erros terms
  matrix[Ng, R] Kgag;
  matrix[Ne, R] Keae;

  //likelihood
  for(i in 1:N){
    y[i] ~ normal(gen[genotype[i]] + env[environment[i]] + xi[genotype[i], environment[i]], sigma[environment[i]]);
  }
  // prioris
  gen ~ normal(Kg*a_g0,sigma_g);
  a_g0 ~ normal(0, sigma_g);
  env ~ normal(0, sigma_e);

  //kernel regression weights
  for(i in 1:Ng){
    for(r in 1:R){
      ag[i,r] ~ normal(0, sigma_ag);
    }
  }

  for(i in 1:Ne){
    for(r in 1:R){
      ae[i,r] ~ normal(0, sigma_ae);
    }
  }


 // latent variables
 Kgag = Kg*ag;
 Keae = Kg*ae;

 for(i in 1:Ng){
   for(r in 1:R){
     hg[i,r] ~ normal(Kgag[i,r], sigma_g);
   }
 }

  for(i in 1:Ne){
   for(r in 1:R){
     he[i,r] ~ normal(Keae[i,r], sigma_e);
   }
 }

 // priors on sigma
 for(j in 1:Ne){
  sigma[j] ~ gamma(alpha,beta);
 }

}
