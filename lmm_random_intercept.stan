data {
  
  int < lower = 1 > N; // Sample size
  vector[N] y; // MIC
  real<lower=0, upper=1> x[N]; //predictor
  int<lower=1> J; //number of PP cluster
  int<lower=1, upper=J> PP[N]; //PP cluster id

}

parameters {
  vector[2] beta; //intercept and slope
  vector[J] u; // varying intercept for PP clusters
  real<lower=0> sigma_e; //error sd
  real<lower=0> sigma_u; //PP cluster sd
}

model {
  real mu;
  //priors

  u ~ normal(0, sigma_u); //PP random effects
  beta[2]~normal(0,0.6);

  // likelihood
  for (i in 1:N){
    
    mu = beta[1] + u[PP[i]] + beta[2] * x[i]; // adds u[cluster[i]] to the mean beta[1](intercept) of distribution of x[i](snp)
    
    y[i] ~ lognormal(mu, sigma_e);
  }
}

