## Run No Interaction Model
library(rstan)
library(bayesplot)

# Function to simulate the data -------------------------------------------

dataSim <- function(Ng, Ne, Kg, Ke, R, sigma_g, sigma_e, sigma_ag, sigma_ae, sigma_y){
  # total of observations
  N <- Ng*Ne
  # simulating genotype
  a_g0 <- rnorm(Ng, 0, sigma_g)
  mu_g <- Kg%*%a_g0
  g <- rnorm(Ng, mu_g, sigma_g)

  # simulating environment
  e <- rnorm(Ne, 0, sigma_e)

  # simulating GxE terms

  #kernel regression weights
  a_g <- matrix(rnorm(Ng*R, 0, sigma_ag), nrow = Ng, ncol = R)
  a_e <- matrix(rnorm(Ne*R, 0, sigma_ae), nrow = Ne, ncol = R)

  # erros terms
  epsilon_Hg <- matrix(rnorm(Ng*R, 0, sigma_g), nrow = Ng, ncol = R)
  epsilon_He <- matrix(rnorm(Ne*R, 0, sigma_e), nrow = Ne, ncol = R)

  # latent variables
  Hg <- Kg%*%a_g + epsilon_Hg
  He <- Ke%*%a_e + epsilon_He
  xi <- Hg%*%t(He)


  x <- expand.grid(1:Ng, 1:Ne)
  names(x) <- c("g", "e")
  x$g <- as.factor(x$g)
  x$e <- as.factor(x$e)

  epsilon_y <- rnorm(Ne, 0, sigma_y)

  y <- g[x[, "g"]] + e[x[, "e"]] + as.vector(xi) + epsilon_y[x[, "e"]]

  y_NI <- g[x[, "g"]] + e[x[, "e"]] + epsilon_y[x[, "e"]]

  return(list(y = y,
              g = g,
              e = e,
              a_g0 = a_g0,
              a_g = a_g,
              a_e = a_e,
              Hg = Hg,
              He = He,
              epsilon_y = epsilon_y,
              x = x,
              y_NI= y_NI,
              Kg = Kg,
              Ke = Ke
  ))


}

Ng = 4
Ne = 6
Kg = matrix(rnorm(Ng*Ng), ncol = Ng)
Ke = matrix(rnorm(Ne*Ne), ncol = Ne)

data <- dataSim(Ng = Ng,
                Ne = Ne,
                Kg = Kg,
                Ke = Ke,
                R = 2,
                sigma_g = 1,
                sigma_e = 1,
                sigma_ag = 1 ,
                sigma_ae = 1,
                sigma_y =  rep(1, Ne))



# Stan code ---------------------------------------------------------------
noIntCode <- "
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
"


# Stan data ---------------------------------------------------------------

noIntData <- list(N = length(data$y_NI),
              Ng = length(data$g),
              Ne = length(data$e),
              genotype = as.integer(data$x$g),
              environment = as.integer(data$x$e),
              y = data$y_NI,
              Kg = data$Kg,
              sigma_g = 1,
              sigma_e = 1,
              alpha = 0.01,
              beta = 0.01
              )


# Stan fit model --------------------------------------------------------------

noIntModel <- stan(model_code = noIntCode,
                   data = noIntData,
                   chains = 1,
                   warmup = 1000,
                   iter = 5000,
                   thin = 1)

print(noIntModel)


# Some plots --------------------------------------------------------------
mcmc_trace(noIntModel, pars = "gen[1]")
mcmc_trace(noIntModel, pars = "env[1]")

posterior <- as.array(noIntModel)
mcmc_areas(
  posterior,
  pars = c("gen[1]", "gen[2]", "gen[3]",
           "env[1]", "env[2]", "env[3]"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)





# Interaction Model -------------------------------------------------------


# Stan code ---------------------------------------------------------------

intCode <- "
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
"

# Stan data ---------------------------------------------------------------

IntData <- list(N = length(data$y),
                Ng = length(data$g),
                Ne = length(data$e),
                genotype = as.integer(data$x$g),
                environment = as.integer(data$x$e),
                y = data$y,
                Kg = data$Kg,
                Kg = data$Ke,
                sigma_g = 1,
                sigma_e = 1,
                alpha = 0.01,
                beta = 0.01
)


# Stan fit model --------------------------------------------------------------

noIntModel <- stan(model_code = noIntCode,
                   data = noIntData,
                   chains = 1,
                   warmup = 1000,
                   iter = 5000,
                   thin = 1)

print(noIntModel)


