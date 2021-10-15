
R <- as.matrix(read.csv("spot-rate(JGB).csv", row.names=1))

T <- nrow(R)  # 時点数
N <- ncol(R)  # 年限数
P <- 4        # VARの次元

# パラメータ
F <- array(0.2,dim=c(P,P))
mu <- array(0,dim=c(P,1))
Sigma_epsilon <- diag(rep(0.1,N))
Sigma_eta <- array(0.1,dim=c(P,P))

# 初期値設定
F<- array(c(0.97,-0.06,-0.01,2.1,-0.04,0.86,-0.09,4.65,0.04,0.003,0.94,0.42,0,0,0,1),dim=c(P,P))
mu<-array(c(0.06,-0.042,-0.02,0.205),dim=c(P,1))
Sigma_epsilon <- diag(rep(0.0001,N))
Sigma_eta <- array(c(0.002,-0.002,-0.002,0,0,0.0005,0.0008,0.025,0,0,0.0034,0.0246,0,0,0,0.08),dim=c(P,P))

z <- array(0,dim=c(T+1,N,4))

# 結果を保存する箱を定義
beta_tt <- array(1,dim=c(P,T+1))
beta_tp1t <- array(1,dim=c(P,T+1))
P_tt <- array(0,dim=c(T+1,P,P))
P_tp1t <- array(0,dim=c(T+1,P,P))
v <- array(0,dim=c(T+1,N))
f <- array(0,dim=c(T+1,N,N))

# ファクターローディングを計算する関数
calcz <- function(beta,N){
  z <- array(0,dim=c(N,4))
  for(m in 1:N){
    z[m,1] <- 1
    z[m,2] <- (1-exp(-beta[4]*m))/(beta[4]*m)
    z[m,3] <- (1-exp(-beta[4]*m))/(beta[4]*m) - exp(-beta[4]*m)
    z[m,4] <- beta[2]*(exp(-beta[4]*m)*(beta[4]*m+1)-1)/(beta[4]*m) + beta[3]*((exp(-beta[4]*m)*(beta[4]*m+1)-1)/(beta[4]*m)+beta[4]*m*exp(-beta[4]*m))
  }
  calcz <- z
}

# 対数尤度関数
calcLikelihood <- function(para){
  calcLikelihood  <- 0
  F <- matrix(para[1:P*P], nrow=P, ncol=P)
  mu <- matrix(para[(P*P+1):(P*(P+1))], nrow=P, ncol=1)
  Sigma_epsilon <- diag(para[(P*(P+1)+1):(P*(P+1)+N)])
  Sigma_eta <- matrix(c(para[P*(P+1)+N+1],para[P*(P+1)+N+2],para[P*(P+1)+N+3],para[P*(P+1)+N+4],0,para[P*(P+1)+N+5],para[P*(P+1)+N+6],
                        para[P*(P+1)+N+7],0,0,para[P*(P+1)+N+8],para[P*(P+1)+N+9],0,0,0,para[P*(P+1)+N+10]),nrow=P, ncol=P)
  
  tryCatch({
    for(t in 2:(T+1)){
      # 予測ステップ
      beta_tp1t[,t] <- (diag(P)-F)%*%mu + F%*%beta_tt[,t-1]
      P_tp1t[t,,] <- F%*%P_tt[t-1,,]%*%F + Sigma_eta

      # 更新ステップ
      z[t,,] <- calcz(beta_tp1t[,t], N)
      v[t,] <- R[t-1,] - z[t,,1:3]%*%beta_tp1t[1:3,t]
      f[t,,] <- z[t,,]%*%P_tp1t[t,,]%*%t(z[t,,]) + Sigma_epsilon
      beta_tt[,t] <- beta_tp1t[,t] + P_tp1t[t,,]%*%t(z[t,,])%*%solve(f[t,,])%*%v[t,]
      P_tt[t,,] <- P_tp1t[t,,]- P_tp1t[t,,]%*%t(z[t,,])%*%solve(f[t,,])%*%z[t,,]%*%P_tp1t[t,,]
      calcLikelihood <- calcLikelihood -N/2*log(2*pi) -1/2*log(det(f[t,,])) -1/2*v[t,] %*%solve(f[t,,])%*%v[t,]
    }
    return(-calcLikelihood)}

    ,error = function(e) {return(100000000000000000)}
    ,warning = function(e) {return(100000000000000000)}
  )
}

para<-c(as.vector(F),as.vector(mu), as.vector(diag(Sigma_epsilon)), c(0.002,-0.002,-0.002,0,0.0005,0.0008,0.025,0.0034,0.0246,0.08))
d<-calcLikelihood(para)

optimizeResult <- optim(par=c(as.vector(F),as.vector(mu), as.vector(diag(Sigma_epsilon)), c(0.002,-0.002,-0.002,0,0.0005,0.0008,0.025,0.0034,0.0246,0.08)), fn=calcLikelihood, method="BFGS")
para <- optimizeResult$par
F <- matrix(para[1:(P*P)], nrow=P, ncol=P)
mu <- matrix(para[(P*P+1):(P*(P+1))], nrow=P, ncol=1)
Sigma_epsilon <- diag(para[(P*(P+1)+1):(P*(P+1)+N)])
Sigma_eta <- matrix(c(para[P*(P+1)+N+1],para[P*(P+1)+N+2],para[P*(P+1)+N+3],para[P*(P+1)+N+4],0,para[P*(P+1)+N+5],para[P*(P+1)+N+6],
                      para[P*(P+1)+N+7],0,0,para[P*(P+1)+N+8],para[P*(P+1)+N+9],0,0,0,para[P*(P+1)+N+10]),nrow=P, ncol=P)

# para<-c(as.vector(F),as.vector(mu), as.vector(diag(Sigma_epsilon)), as.vector(Sigma_eta))
# d<-calcLikelihood(para)
# optimizeResult <- optim(par=c(as.vector(F),as.vector(mu), as.vector(diag(Sigma_epsilon)), as.vector(Sigma_eta)), fn=calcLikelihood, method="BFGS")


# # 対数尤度関数
# calcLikelihood <- function(param1, param2){
#   fpara <- array(0,dim=c(T,N,N))
#   vpara <- array(0,dim=c(T,N))
#   for(i in 1:T){
#     for(j in 1:N){
#       fpara[i,j,] <- param1[(N*N(i-1)+N*(j-1)+1):(N*N(i-1)+N*(j-1)+N)]
#     }
#     vpara[i,] <- param2[(N*(i-1)+1):(N*(i-1)+N)]
#     calcLikelihood <- calcLikelihood -N/2*log(2*pi) -1/2*log(det(fpara[i,,])) -1/2*vpara[i,] %*%solve(fpara[i,,])%*%t(vpara[i,])
#   }
# }