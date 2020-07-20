data {
  int T;             // データ数
  vector[4] pi[T];   // 観測データ
}

parameters{
  vector<lower=-1,upper=1>[4] rho_tau;
  matrix[4,4] b1;
  matrix[4,4] b2;
  vector<lower=0>[4] sigma_z;
  real<lower=0> sigma_tau;
  vector[T] tau;
  vector[4] z[T];
}

transformed parameters{
  cov_matrix[4] Sigma_z;    // zの分散共分散行列（対角行列）
  vector[4] tau_vector[T];
  vector[4] Sigma_tau;
  cov_matrix[4] Sigma_pi;
  
  Sigma_tau[1] = sigma_tau;
  Sigma_tau[2] = sigma_tau;
  Sigma_tau[3] = sigma_tau;
  Sigma_tau[4] = sigma_tau;
  Sigma_z[1,1] = sqrt(sigma_z[1]^2 + sigma_tau^2*rho_tau[1]^2);
  Sigma_z[2,2] = sqrt(sigma_z[2]^2 + sigma_tau^2*rho_tau[2]^2);
  Sigma_z[3,3] = sqrt(sigma_z[3]^2 + sigma_tau^2*rho_tau[3]^2);
  Sigma_z[4,4] = sqrt(sigma_z[4]^2 + sigma_tau^2*rho_tau[4]^2);
  Sigma_z[1,2] = 0;
  Sigma_z[1,3] = 0;
  Sigma_z[1,4] = 0;
  Sigma_z[2,1] = 0;
  Sigma_z[2,3] = 0;
  Sigma_z[2,4] = 0;
  Sigma_z[3,1] = 0;
  Sigma_z[3,2] = 0;
  Sigma_z[3,4] = 0;
  Sigma_z[4,1] = 0;
  Sigma_z[4,2] = 0;
  Sigma_z[4,3] = 0;
  Sigma_pi[1,1] = sqrt(sigma_tau^2+sigma_z[1]^2+sigma_tau^2*rho_tau[1]^2);
  Sigma_pi[2,2] = sqrt(sigma_tau^2+sigma_z[2]^2 + sigma_tau^2*rho_tau[2]^2);
  Sigma_pi[3,3] = sqrt(sigma_tau^2+sigma_z[3]^2 + sigma_tau^2*rho_tau[3]^2);
  Sigma_pi[4,4] = sqrt(sigma_tau^2+sigma_z[4]^2 + sigma_tau^2*rho_tau[4]^2);
  Sigma_pi[1,2] = 0;
  Sigma_pi[1,3] = 0;
  Sigma_pi[1,4] = 0;
  Sigma_pi[2,1] = 0;
  Sigma_pi[2,3] = 0;
  Sigma_pi[2,4] = 0;
  Sigma_pi[3,1] = 0;
  Sigma_pi[3,2] = 0;
  Sigma_pi[3,4] = 0;
  Sigma_pi[4,1] = 0;
  Sigma_pi[4,2] = 0;
  Sigma_pi[4,3] = 0;

  for(i in 1:T){
    tau_vector[i] = rep_vector(tau[i],4);
  }
}

model{
  // 事前分布の設定
  rho_tau ~ multi_normal(rep_vector(0,4),diag_matrix(rep_vector(40,4)));
  to_vector(b1) ~ multi_normal(rep_vector(0.5,16),diag_matrix(rep_vector(0.02,16)));
  to_vector(b2) ~ multi_normal(rep_vector(0,16),diag_matrix(rep_vector(0.02,16)));
  sigma_z[1] ~ inv_gamma(5,2.5);
  sigma_z[2] ~ inv_gamma(5,2.5);
  sigma_z[3] ~ inv_gamma(5,2.5);
  sigma_z[4] ~ inv_gamma(5,2.5);
  sigma_tau ~ inv_gamma(2.5,1);
  
  for(i in 3:T){
    tau[i] ~ normal(tau[i-1], sigma_tau);
    z[i] ~ multi_normal(b1*z[i-1]+b2*z[i-2], Sigma_z);
    pi[i] ~ multi_normal(tau_vector[i-1]+b1*z[i-1]+b2*z[i-2],Sigma_pi);
  }
}



