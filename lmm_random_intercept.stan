data {
  real<lower=0, upper=1> h_squared; //heritability
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
}

transformed_parameters {
  real<lower=0> sigma_u; //PP cluster sd
  sigma_u = sigma_e * 1/sqrt(h_squared - 1);
}

model {
  real mu;
  //priors

  u ~ normal(0, sigma_u); //PP random effects

  // likelihood
  for (i in 1:N){

    mu = beta[1] + u[PP[i]] + beta[2] * x[i]; // adds u[cluster[i]] to the mean beta[1](intercept) of distribution of x[i](snp)

    y[i] ~ normal(mu, sigma_e);
  }
}

