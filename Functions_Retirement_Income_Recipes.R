
# 
RGOA <- function(g,v,N){
  # 成長率gの年金における現在価値係数を計算する関数
  #
  # g：成長率
  # v：CFの割引率
  # N 期間数
  (1-((1+g)^N)*(1+v)^(-N))/((v-g)/(1+g))
}

SMCR<-function(x,f,w,g,v,R,D){
  # 最適な消費額（一定）を計算する関数
  #
  # x：現在の年齢
  # f：現在の金融資本
  # w：現在の賃金
  # g：賃金の成長率（賃金が年あたりなら年次）
  # v：CFの評価率（賃金が年あたりなら年次）
  # R：定年年齢
  # D：死亡年齢
  (w*RGOA(g,v,R-x)+f)/RGOA(0,v,D-x)
}

OLCF<-function(x,x0,f0,w0,g,v,R,D){
  # 特定の年齢x歳で持つべき最適な金融資本を計算する関数
  #
  # x：現在の年齢
  # x0：初期年齢（働き始めor貯蓄始めた年齢）x0<=x
  # f0：x0での金融資本
  # w0：w0での賃金
  # g：賃金の成長率（賃金が年あたりなら年次）
  # v：CFの評価率（賃金が年あたりなら年次）
  # R：定年年齢
  # D：死亡年齢
  if (x<=R){
    SMCR(x0,f0,w0,g,v,R,D)*RGOA(0,v,D-x)-
      (w0*(1+g)^(x-x0))*RGOA(g,v,R-x)
  }else {
    SMCR(x0,f0,w0,g,v,R,D)*RGOA(0,v,D-x)
  }
}

FCMW<-function(x,x0,g,v,R,D){
  # 年齢x歳で持つべき金融資本の賃金に対する比率を計算する関数
  #
  # x：現在の年齢
  # x0：初期年齢（働き始めor貯蓄始めた年齢）x0<=x
  # g：賃金の成長率（賃金が年あたりなら年次）
  # v：CFの評価率（CFが年あたりなら年次）
  # R：定年年齢
  # D：死亡年齢
  ((1+g)^(x0-x))*RGOA(g,v,R-x0)*RGOA(0,v,D-x)/RGOA(0,v,D-x0)-RGOA(g,v,R-x)
}

PL<-function(v,c,F){
  # ポートフォリオの寿命を計算する関数
  # 
  # v：評価率
  # c：消費額
  # F：金融資本
  if (c/F <= v){
    999　# 無限大の代わり
  }else{
    if (v==0){
      F/c
    }
    else{
      (1/v)*log( (c/F)/(c/F-v) )
    }
  }
}

DTRJ<-function(t,v,c,F){
  # 任意の時点におけるの金融資本を計算する関数
  # 
  # t：時点
  # v：評価率
  # c：消費額
  # F：金融資本
  F*exp(v*t)-c*(exp(v*t)-1)/v
}

PLSM<-function(F,c,nu,sigma,N){
  # ポートフォリオの寿命をシミュレーションする関数  
  #
  # F：初期の金融資本
  # c：消費額
  # nu：投資収益率の平均
  # sigma：投資収益率の標準偏差
  # N：シミュレーション回数
  path<-matrix(nrow=N,ncol=100)
  PLV<-c()
  for (i in 1:N){
    return<-exp(rnorm(100,nu,sigma))
    path[i,1]<-F
    for (j in 2:100){
      path[i,j]<-path[i,j-1]*return[j]-c
      if (path[i,j]<=0) {break}
    }
    PLV[i]=j
  }
  PLV
}

PLSM.SR<-function(F,c,nu,sigma,N){
  # 10年間の幾何平均を保存し，ポートフォリオの寿命をシミュレーションする関数
  #
  # F：初期の金融資本
  # c：消費額
  # nu：投資収益率の平均
  # sigma：投資収益率の標準偏差
  # N：シミュレーション回数
  path<-matrix(nrow=N,ncol=100)
  PLM<-matrix(nrow=N,ncol=4)
  for (i in 1:N){
    return<-exp(rnorm(100,nu,sigma))
    PLM[i,1]<-prod(return[1:10])^(1/10)-1
    PLM[i,2]<-prod(return[11:20])^(1/10)-1
    PLM[i,3]<-prod(return[21:30])^(1/10)-1
    path[i,1]<-F
    for (j in 2:100){
      path[i,j]<-path[i,j-1]*return[j]-c
      if (path[i,j]<=0) {break}
    }
    PLM[i,4]=j
  }
  PLM
}

SPQR<-function(x,y,qx){
  # 1年死亡率ベクトルからn年生存率を計算する関数
  # 
  # x：現在の年齢（条件）
  # y：求めたい生存率の年齢
  # 死亡率ベクトル
  LT<-cumprod(1-qx)
  LT<-c(1,LT)
  LT[y+1]/LT[x+1]
}

TPXG<-function(x,t,m,b){
  # ゴンぺルツ仮定の下，生存率を求める関数
  # 
  # x：現在の年齢
  # t：時間
  # m：ゴンぺルツパラメータ，最頻値
  # b：ゴンぺルツパラメータ，分散係数
  exp(exp((x-m)/b)*(1-exp(t/b)))
}

lambda<-function(t){
  # 数値積分するための，ゴンぺルツ仮定の下，生存率を求める関数
  # 
  # t：時間
  -(1/b)*exp((x+t-m)/b)
}

GRAN<-function(N,x,m,b){
  # ゴンぺルツ乱数（時点のシミュレーション）を生成する関数
  # 
  # N：シミュレーション回数
  # x：現在の年齢
  # m：ゴンぺルツパラメータ，最頻値
  # b：ゴンぺルツパラメータ，分散係数
  b*log(1-log(runif(N))*exp((m-x)/b))
}

LTLD<-function(z){
  # Gompertzr乱数によりシミュレーションされた寿命からコホート生命表を作成する関数
  # 
  # z：死亡時点ベクトル
  z<-sort(z); N<-length(z); LT<-c();
  LT[1]<-N; T<-ceiling(max(z))+1
  for (i in 2:T){
    LT[i]<-N-length(z[z<=i-1])
  }
  LT
}

#Remember to set x,b and qx values.
gap<-function(m){
  # 最頻値を与えて，ゴンペルツ分布の解析解とシミュレーションによる数値解の絶対偏差の和を計算する関数
  # 
  # m：最頻値
  gap<-c()
  for (i in 1:45){
    gap[i]<-abs(TPXG(x,i,m,b)-SPQR(x,x+i,qx))}
  sum(gap)
}

LRPG<-function(v,xi,x,m,b){
  # 
  #
  # v：CFの評価率（CFが年あたりなら年次）
  # xi：支出率
  # x：現在の年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  if (xi <= v){
    ruin<-0
  }else{
    ruin<-exp( exp((x-m)/b) * (1 - (xi/(xi-v))^(1/(v*b))))
  }
  ruin
}

VARPHI.SM<-function(N,x,m,b,xi,nu,sigma){
  # LRPをシミュレーションするための関数(9.7節）
  # 
  # N：シミュレーション回数
  # x：現在の年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # xi：支出率
  # nu：投資収益率の期待値
  # sigma：投資収益率の標準偏差
  V<-c()
  wks<-round(GRAN(N,x,m,b)*52,0)+1
  for (i in 1:N){
    sB<-sigma*sqrt(1/52)*cumsum(rnorm(wks[i]))
    Z<-exp((nu/52)*c(1:(wks[i]))+sB)
    V[i]<-sum((1/Z)/52)
  }
  sum(V>=1/xi)/length(V)
}

G<-function(a,c){
  # 
  # 
  # a
  # c
  #Avoid c=0, when a=0.
  integrand<-function(t){
    (t^(a-1)*exp(-t))
  }
  integrate(integrand,c,Inf)$value
}

a<-function(v,x,m,b){
  # 
  # 
  # v：CFの評価率（CFが年あたりなら年次）
  # x：現在の年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  b*exp(exp((x-m)/b)+(x-m)*v)*G(-b*v,exp((x-m)/b))
}

VARPHI.MM<-function(x,m,b,xi,nu,sigma){
  # 
  # 
  # x：現在の年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # xi：支出率
  # nu：投資収益率の期待値
  # sigma：投資収益率の標準偏差
  mu<-nu+(0.5)*sigma^2
  M1<-a(mu-sigma^2,x,m,b)
  M2<-(a(mu-sigma^2,x,m,b)
       -a(2*mu-3*sigma^2,x,m,b))/(mu/2-sigma^2)
  alpha<-(2*M2-M1^2)/(M2-M1^2)
  beta<-(M2-M1^2)/(M2*M1)
  #Shape is Alpha, Scale is Beta
  pgamma(xi,shape=alpha,scale=beta,lower.tail=TRUE)
}

GILA<-function(x,v,m,b){
  # 即時生命年金の価値を評価する
  #
  # x：現在の年齢
  # v：評価率
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  omega<-x+b*log(1+10*log(10)*exp((m-x)/b));      # 最大生存年齢
  dt<-1/52; grid<-ceiling((omega-x)/dt); t<-(1:grid)*dt    # grid：グリッド数，t：週次の時点
  pgrid<-exp(exp((x-m)/b)*(1-exp(t/b)))  # 生存率（8.7式）
  rgrid<-exp(-v*t)
  sum(pgrid*rgrid)*dt
}

GTLA<-function(x,tau,v,m,b){
  # 一時年金の価値を評価する
  #
  # x：現在の年齢
  # tau：年金期間
  # v：評価率
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  if (tau==0){0}
  else {
    dt<-1/52; grid<-ceiling(tau/dt); t<-(1:grid)*dt
    pgrid<-exp(exp((x-m)/b)*(1-exp(t/b)))
    rgrid<-exp(-v*t)
    sum(pgrid*rgrid)*dt}
}

GDLA<-function(x,y,v,m,b){
  # 据置年金の価値を評価する
  #
  # x：現在の年齢
  # y：年金開始年齢
  # v：評価率
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  GILA(x,v,m,b)-GTLA(x,(y-x),v,m,b)
}

IDDR<-function(Fx,pi,x,m,b,r,rho,gam){
  # インテリジェントドローダウンを計算する関数（c0-π)/Fx 
  #
  # Fx：金融資本
  # pi：年金収入
  # x:年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # r：金利
  # rho：主観的割引率
  # gam：長寿リスク回避度
  # Positive Pension (pi>0) Only.
  k<-(r-rho)/gam
  WDT<-function(tau){
    # 富が枯渇する時点を計算する関数（c*t/π-1を計算）
    K1<-((Fx/pi+1/r)*exp(r*tau)-1/r)/
      (GTLA(x-b*log(gam),tau,r-k,m,b)*exp(r*tau))
    K2<-exp(k*tau)
    K3<-TPXG(x,tau,m,b)^(1/gam)
    K1*K2*K3-1
  }
  tau<-uniroot(WDT,lower=0,upper=100)$root
  Cx<-((Fx+pi/r)*exp(r*tau)-pi/r)/
    (GTLA(x-b*log(gam),tau,r-k,m,b)*exp(r*tau))
  (Cx-pi)/Fx
}

IDDR0<-function(Fx,x,m,b,r,rho,gam){
  # 年金収入がない場合のインテリジェントドローダウンを計算する関数（約分できてシンプルに）
  # 
  # Fx：金融資本
  # x:年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # r：金利
  # rho：主観的割引率
  # gam：長寿リスク回避度
  k<-(r-rho)/gam
  Cx<-Fx/GILA(x-b*log(gam),r-k,m,b)
  Cx/Fx
}

WDT.PSI<-function(psi,x,m,b,r,rho,gam){
  # 富の枯渇時間を求めるための関数
  # 
  # psi：個人のバランスシートの内，年金化された資産の割合
  # x:年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # r：金利
  # rho：主観的割引率
  # gam：長寿リスク回避度
  f<-function(tau){
    k<-(r-rho)/gam
    K1<-(((1-psi)*GILA(x,r,m,b)/psi+1/r)*exp(r*tau)-1/r)/
      (GTLA(x-b*log(gam),tau,r-k,m,b)*exp(r*tau))
    K2<-exp(k*tau)
    K3<-TPXG(x,tau,m,b)^(1/gam)
    K1*K2*K3-1
  }
  uniroot(f,lower=0,upper=100)$root
}

PSI.OPT<-function(tau,x,m,b,r,rho,gam){
  # WDTをインプットとして，psiを求める関数
  # 
  # tau：WDT
  # x:年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # r：金利
  # rho：主観的割引率
  # gam：長寿リスク回避度
  k<-(r-rho)/gam
  a1<-GILA(x,r,m,b)
  a2<-GTLA(x-b*log(gam),tau,r-k,m,b)
  K1<-(a2/a1)*exp(-k*tau)*TPXG(x,tau,m,b)^(-1/gam)
  K2<-(1/r)*(1/a1)*(exp(-r*tau)-1)
  1/(K1+K2+1)
}

af<-function(v,x,m,b){
  # 
  # 
  # v
  # x:年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  b*exp(exp((x-m)/b)+(x-m)*v)*G(-b*v,exp((x-m)/b))
}
# Converts Utility Benefit of Pension Annuity into Dollars.
DLTA<-function(x,m,b,v,gam){
  # 
  # x:年齢
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  # v
  # gam
  (af(v,x,m,b)/af(v,x-b*log(gam),m,b))^(gam/(1-gam))-1
}

ANHG<-function(x,v,h0,g){
  # The Gompertz Annuity Valuation Model (GAVM)
  # Using equation (A.33) from Milevsky (JPEF, 2020)
  hx<-h0*exp(g*x)
  G(-v/g,hx/g)/(g*exp((-1/g)*(hx+v*log(hx/g))))
}
UDHG<-function(x,h0,g,v,gam){
  # The Delta from Annuitization
  # This is based on equation (12.13)
  (ANHG(x,v,h0,g)/ANHG(x-log(gam)/g,v,h0,g))^(gam/(1-gam))-1 
}

#Gompertz Makeham Survival Probability
GMSP<-function(x,t,lam,m,b){
  # ゴンペルツモデルにおける生存確率を計算する関数
  # 
  # x:年齢
  # t：時点
  # lam：ハザード率
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  exp(-lam*t+exp((x-m)/b)*(1-exp(t/b)))
}

# The Gompertz Makeham Probability Density Function
GMPDF<-function(x,t,lam,m,b){
  # ゴンペルツ・マケハムモデルの確率密度関数を計算する関数
  # 
  # x:年齢
  # t：時点
  # lam：マケハム定数
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  exp(-lam*t+exp((x-m)/b)*(1-exp(t/b)))*(lam+(1/b)*exp((x+t-m)/b))
}

# Gompertz Makeham Life Expectancy
GMLE<-function(x,lam,m,b){
  # ゴンペルツ・マケハムモデルにおける余命の期待値を計算する関数
  # 
  # x:年齢
  # lam：マケハム定数
  # m：ゴンペルツパラメータ（最頻値）
  # b：ゴンペルツパラメータ（分散係数）
  b*G(-lam*b,exp((x-m)/b))/exp((m-x)*lam-exp((x-m)/b))
}

CLAM<-function(b,x_star,lam_star){
  # ゴンペルツ・マケハムモデルの下でのmの近似値を計算する関数
  # 
  # b：ゴンペルツパラメータ（分散係数）
  # x_star：年齢
  # lam_star：マケハム定数
  x_star-b*(log(b*lam_star))
}

LRAG<-function(x,x_star,g_hat,g,lam_hat,lam,lam_star){
  # リスク調整後年齢を計算する関数
  #
  # x：現在の年齢
  # x_star：ハザードレートが一定となる年齢
  # g_hay：gのグローバル平均
  # g：国固有のg
  # lam_hat：λのグローバル平均
  # lam：国固有のλ
  # lam_star：ハザードレートが一定となった後の値
  x_star+(1/g_hat)*
    log(exp(g*(x-x_star))-(lam_hat-lam)/lam_star)
}