library(rstan)
library(zoo)
library(xts)
library(ggfortify)
library(ggplot2)
library(bayesplot)

# 計算の高速化
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

# 消費者物価指数データ読み込み
CPI <- read.csv("CPI.csv",header=TRUE)
CPI <- as.xts(read.zoo(CPI))

start_date <- "1983-03-01"
end_date <- "2016-03-01"

pi <- CPI['1983-03-01::2016-03-01',1]

test <- read.csv("test.csv",header=TRUE)
T <- nrow(test)

# データの準備
data_list <- list(
  T = T,
  pi = test$生鮮食品を除く総合
)

# モデルの推定
BN_stan <- stan(
  file = "BN.stan",
  data = data_list,
  seed = 1111,
  iter = 20000,
  warmup = 10000,
  chains = 4,
  control = list(adapt_delta=0.9,max_treedepth=15)
)

# 推定結果の確認
print(BN_stan,probs=c(0.25,0.5,0.975),
      pars=c("b1","b2","sigma_z","sigma_tau","rho_tau"))
traceplot(BN_stan,pars=c("b1","b2"))#,inc_warmup=T)
traceplot(BN_stan,pars=c("sigma_z","sigma_tau"))
traceplot(BN_stan,pars=c("rho_tau"))

stan_dens(BN_stan,pars="b1",separate_chains = T)
stan_dens(BN_stan,pars="b2",separate_chains = T)
stan_dens(BN_stan,pars="rho_tau",separate_chains = T)
stan_dens(BN_stan,pars="sigma_z",separate_chains = T)
stan_dens(BN_stan,pars="sigma_tau",separate_chains = T)

# MCMCサンプルの抽出
mcmc_sample <- rstan::extract(BN_stan,permuted = FALSE)

b1 <- mcmc_sample[,,"b1"]
b2 <- mcmc_sample[,,"b2"]
rho_tau <- mcmc_sample[,,"rho_tau"]
sigma_z <- mcmc_sample[,,"sigma_z"]
sigma_tau <- mcmc_sample[,,"sigma_tau"]

autoplot(ts(b1),facets=F,#4つのセットを同時にグラフ
         ylib="b1",main="トレースプロット")

# ggplotによる図示  
b1_df <- data.frame(b1 = as.vector(b1))
ggplot(data=b1_df, mapping=aes(x=b1))+geom_density(size=1.5)

sigma_tau_df <- data.frame(sigma_tau = as.vector(sigma_tau))
ggplot(data=sigma_tau_df, mapping=aes(x=sigma_tau))+geom_density(size=1.5)

# bayesplotによる事後分布の可視化
mcmc_sample_para <- list(b1,b2)
mcmc_hist(mcmc_sample_para,pars=c("b1","b2"))
mcmc_combo(mcmc_sample_para)
