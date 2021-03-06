# cross validation
# population minimizer
# inverstigation of asymptotic variances
rm (list =ls()); gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
library(e1071)
library(quadprog)
library(MASS)
library(igraph)
library(glmnet)
library(dplyr)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
f = function(x,y,df)  pt(x-y, df = df)*dt(y,df = df)


# set parameters in simulation
nvec = seq(25,500,by = 25)
gamma.vec = seq(0,1,by = 0.25)
df = 1 ;     alpha = 5
lambda.vec = c(1.2,1.1,0.6,0)*4
#lambda.vec = c(1.2,0.8,0.6,0)*4
pop_list = list()
pop_beta_mat = NULL
rankcons.mat1 = rankcons.mat2 = rankcons.mat3 = NULL
cov.list = list()
p = length(lambda.vec)
Gmat = Qmat = matrix(0,p,p)
kn = 3
k_fold = 5
for (i in 1:(p-1) )
{
  for (j in (i+1):p)
  {
    un = kn
    if ( (i == 2 & j == 3) ) un = alpha*kn
    if ( (i == 1 & j == 4) ) un = alpha*kn
    if ((i==3) & (j==4)) un = alpha*kn
    tmp.vec = rep(0,p)
    tmp.vec[i] = 1 ; tmp.vec[j] = -1
    #v = pt((lambda.vec[i]- lambda.vec[j])/2, df = df)
    v = integrate(f, -Inf, Inf, subdivisions = 1e+5L,
                  x = (lambda.vec[i] - lambda.vec[j]), df = df,
                  abs.tol = 1e-8)$value
    #v = plogis((lambda.vec[i]- lambda.vec[j]))
    Gmat[i,j] = v 
    Gmat[j,i] = 1-v 
    Qmat[i,j] = Qmat[j,i] = un
  }
}
#### complete the set-up here!



a1 = c()
a2 = c()
a3 = c()
for ( ii in 1:500)
{
  set.seed(ii)
  Gmat_obs = gen.Gmat_obs(Gmat, Qmat)
  tn = sum(Qmat)/2
  cv_m = cv_mat_fun(Gmat_obs, Qmat, k_fold = k_fold) 
  
  cv.vec = rep(0, length(gamma.vec))
  for (u in 1:length(gamma.vec))
  {
    gamma.v = gamma.vec[u]
    vv = 0
    for (k in 1:k_fold)
    {
      cv_tr = cv_m[cv_m[,4]!=k,]
      cv_te = cv_m[cv_m[,4]==k,]
      tb_tr = cv_table_fun(cv_tr)
      tb_te = cv_table_fun(cv_te)
      Gh  = tb_tr$G/tb_tr$Q
      Gh[is.na(Gh) ] = 0
      Q = tb_tr$Q
      # set up Qpmat 
      Qp = Q/sum(Q)*2
      # set up (x,y) : complete comparisons
      dfit = gen.designR(p)
      Qmat_fit = list()
      Qmat_fit$Qpmat = Qp
      Qmat_fit$Gmat_hat = Gh
      Qmat_fit$x = dfit$x
      Qmat_fit$y = dfit$y
      # minimizer
      bt_data = Qmat_fit
      Qp = bt_data$Qpmat
      wmat = (1/Qp)^gamma.v ; diag(wmat) = 0
      wmat[wmat==Inf] = 0
      bt_data$Qpmat = Qp*wmat
      try_fit = try({beta_tilde = btFun(bt_data)}, silent = T)
      if (class(try_fit) == 'try-error') next
      
      # test set and compute pred err
      Gmat_te  = tb_te$G/tb_te$Q 
      Gmat_te[is.na(Gmat_te) ] = NA
      
      v = 0
      for (i in 1:(p-1))
      {
        for (j in (i+1):p)
        {
          if (is.na(Gmat_te[i,j]))
          {
            next
            v = v + 0.5
          }
          
          if (Gmat_te[i,j]==0.5)
          {
            next
            v = v + 0.5
          }
          
          if (sign((beta_tilde[i] - beta_tilde[j])*(Gmat_te[i,j]-0.5))>0)
          {
            v = v + 1
          }
          
        }
      }
      vv = vv + v
    }
    cv.vec[u] = vv
  }
  
  opt_idx= max(which(cv.vec == max(cv.vec)))
  
  Gmat_hat = Gmat_obs/Qmat ; diag(Gmat_hat) = 0
  Qmat_hat = Qmat/sum(Qmat)*2
  Qmat_fit = list()
  Qmat_fit$Qpmat = Qmat_hat
  Qmat_fit$Gmat_hat = Gmat_hat
  Qmat_fit$x = dfit$x
  Qmat_fit$y = dfit$y
  # cv minimizer
  gamma.v = gamma.vec[opt_idx]
  bt_data = Qmat_fit
  Qp = bt_data$Qpmat
  wmat = (1/Qp)^gamma.v ; diag(wmat) = 0
  wmat[wmat==Inf] = 0
  bt_data$Qpmat = Qp*wmat

  try_fit = try({beta_tilde = btFun(bt_data)}, silent = T)
  if (class(try_fit) == 'try-error') 
  {
    next
  }
  a1[ii] = sum(order(beta_tilde) == 4:1) == 4
  
  
  # BT 
  gamma.v = 0
  bt_data = Qmat_fit
  Qp = bt_data$Qpmat
  wmat = (1/Qp)^gamma.v ; diag(wmat) = 0
  wmat[wmat==Inf] = 0
  bt_data$Qpmat = Qp*wmat
  beta_tilde = btFun(bt_data)
  a2[ii] = sum(order(beta_tilde) == 4:1) == 4
  
  
  # gBT 
  gamma.v = 1
  bt_data = Qmat_fit
  Qp = bt_data$Qpmat
  wmat = (1/Qp)^gamma.v ; diag(wmat) = 0
  wmat[wmat==Inf] = 0
  bt_data$Qpmat = Qp*wmat
  beta_tilde = btFun(bt_data)
  a3[ii] = sum(order(beta_tilde) == 4:1) == 4
  cat(ii,'\n')
}
# iteration
mean(a1) ; mean(a2) ; mean(a3)

