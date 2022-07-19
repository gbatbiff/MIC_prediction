data {
  int < lower = 1 > N; // Sample size
  int<lower=0, upper=1> x[N]; //predictor
  int<lower=1> J; //number of PP cluster
  int<lower=1, upper=J> PP[N]; //PP cluster id
  vector[N] y; // MIC
  
}

parameters {
  vector[2] beta; //intercept and slope
  vector[J] u; // varying intercept for PP clusters
  real<lower=0> sigma_e; //error sd
  real<lower=0> sigma_u; //PP cluster sd
}


model {
  vector[N] mu_PP;
  //priors
  beta[2]~normal(0,0.6); //sigma=h²
  u ~ normal(0, sigma_u); //PP clusters random effects
     
  // likelihood 
   for(i in 1:N){
      
      mu_PP[i] = beta[1] + u[PP[i]] + beta[2] * x[i]; 
 	    y[i] ~ lognormal(mu_PP[i], sigma_e);
 	  
 	 
  }
}
