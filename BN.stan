data {
  int T;
  vector[T] pi;
}

parameters{
  real rho_tau;
  real b1;
  real b2;
  real<lower=0> sigma_z;
  real<lower=0> sigma_tau;
  vector[T] tau;
  vector[T] z;
}

model{
  // –‘O•ª•z‚Ìİ’è
  rho_tau ~ normal(0,40);
  b1 ~ normal(0.5,0.02);
  b2 ~ normal(0,0.02);
  sigma_z ~ inv_gamma(5,2.5);
  sigma_tau ~ inv_gamma(2.5,1);
  
  tau[2:T] ~ normal(tau[1:T-1], sigma_tau);
  z[3:T] ~ normal(z[2:T-1]*b1+z[1:T-2]*b2, rho_tau+sigma_z);
  pi ~ normal(tau+z,sigma_tau+sigma_z);
}



