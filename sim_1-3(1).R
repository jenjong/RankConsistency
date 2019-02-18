# population minimizer
# inverstigation of asymptotic variances
rm (list =ls()); gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
library(e1071)
library(quadprog)
library(MASS)
library(igraph)
library(glmnet)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
f = function(x,y,df)  pt(x-y, df = df)*dt(y,df = df)


# set parameters in simulation
nvec = seq(25,500,by = 5)
gamma.vec = seq(0,1,by = 0.01)
df = 1;     alpha = 10
#lambda.vec = c(1.2,1.1,0.6,0)*4
lambda.vec = c(1.2,0.8,0.6,0)*4
pop_list = list()
pop_beta_mat = NULL
rankcons.mat1 = rankcons.mat2 = rankcons.mat3 = NULL
cov.list = list()
for (iter in 1:length(gamma.vec))
{
  gamma.v = gamma.vec[iter]
  pop_min = list()
  {
    p = length(lambda.vec)
    idx = 1
    
    Gmat = Qmat = matrix(0,p,p)
    for (i in 1:(p-1) )
    {
      for (j in (i+1):p)
      {
        un = 1
        if ( (i == 2 & j == 3) ) un = alpha
        if ( (i == 1 & j == 4) ) un = alpha
        if ((i==3) & (j==4)) un = alpha
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
        idx = idx + 1
      }
    }
    # set up Qpmat 
    Qpmat = Qmat/sum(Qmat)*2
    
    # set up (x,y) : complete comparisons
    x = matrix(0, p*(p-1), p)
    y = rep(0, p*(p-1) )
    ix = 1
    for (i in 1:p)
    {
      for (j in 1:p)
      {
        if (i == j) next
        jx1 = min(i,j)
        jx2 = max(i,j)
        x[ix,jx1] = 1; x[ix,jx2] = -1
        if (i<j) y[ix] = 1
        ix = ix + 1
      }
    }
    x = x[,-p]
    
    Qmat_fit = list()
    Qmat_fit$Qpmat = Qpmat
    Qmat_fit$Gmat_hat = Gmat
    Qmat_fit$x = x
    Qmat_fit$y = y
    
    # population minimizer
    bt_data = Qmat_fit
    Qpmat = bt_data$Qpmat
    wmat = (1/Qpmat)^gamma.v ; diag(wmat) = 0
    bt_data$Qpmat = Qpmat*wmat
    beta_tilde = btFun(bt_data)
    pop_beta_mat = rbind(pop_beta_mat, beta_tilde)
    pop_min$est = beta_tilde 
    
    
    # true probability
    true_prob = Gmat
    # variance of score function
    v_s = matrix(0,p-1,p-1)
    # diag
    for ( i in 1:(p-1))
    {
      v_s[i,i] = sum(true_prob[i,-i]*(1-true_prob[i,-i])*Qpmat[i,-i]*
                       wmat[i,-i]^2)
    }
    # off-diag
    for ( i in 1:(p-1))
    {
      for (j in 1:(p-1))
      {
        if (i==j) next
        v_s[i,j] = - true_prob[i,j]*(1-true_prob[i,j])*Qpmat[i,j]*
          wmat[i,j]^2
      }
    }
    
    # estimated probability
    est_prob = matrix(0,p,p)
    for (i in 1:p)
    {
      for (j in 1:p)
      {
        if (i==j) next
        v = exp(beta_tilde[i]-beta_tilde[j])
        est_prob[i,j] = v/(1+v)
      }
    }
    
    # Hessian
    H = matrix(0,p-1,p-1)
    # diag
    for ( i in 1:(p-1))
    {
      H[i,i] = sum(est_prob[i,-i]*(1-est_prob[i,-i])*Qpmat[i,-i]*
                     wmat[i,-i])
    }
    # off-diag
    for ( i in 1:(p-1))
    {
      for (j in 1:(p-1))
      {
        if (i==j) next
        H[i,j] = - est_prob[i,j]*(1-est_prob[i,j])*Qpmat[i,j]*
          wmat[i,j]
      }
    }
    # Asymptotic variance 
    FH = solve(H)%*%v_s%*%solve(H)
    pop_min$cov = FH
    cov.list[[iter]] = FH
  }
  pop_list[[iter]] = pop_min
  FH = pop_min$cov
  beta_tilde = pop_min$est
  #
  s1 = FH[1,1] + FH[2,2] - 2*FH[1,2]
  s2 = FH[2,2] + FH[3,3] - 2*FH[2,3]  
  s3 = FH[3,3]
  m1= beta_tilde[1] - beta_tilde[2]
  m2= beta_tilde[2] - beta_tilde[3]
  m3= beta_tilde[3]
  
  rankcons.prob = c()
  for (j in 1:length(nvec))
  {
    rankcons.prob[j] = 1- pnorm(0, m1, sqrt(s1)/sqrt(nvec[j]))
  }
  rankcons.mat1 = rbind(rankcons.mat1, rankcons.prob)
  
  
  rankcons.prob = c()
  for (j in 1:length(nvec))
  {
    rankcons.prob[j] = 1- pnorm(0, m2, sqrt(s2)/sqrt(nvec[j]))
  }
  rankcons.mat2 = rbind(rankcons.mat2, rankcons.prob)
  
  
  rankcons.prob = c()
  for (j in 1:length(nvec))
  {
    rankcons.prob[j] = 1- pnorm(0, m3, sqrt(s3)/sqrt(nvec[j]))
  }
  rankcons.mat3 = rbind(rankcons.mat3, rankcons.prob)
  
}

# figure for Prob(lambda_j > lambda_k)
  rankcons.mat = rankcons.mat3
  rownames(rankcons.mat) = NULL
  colnames(rankcons.mat) = nvec
  rownames(rankcons.mat) = gamma.vec
  library(ggplot2)
  library(reshape)
  rmat = rankcons.mat
  #lim.v = 0.3
  #rmat[rmat<=lim.v]  = lim.v
  fig_data = melt(rmat)
  fig <- ggplot(fig_data , aes(x = X2, y = X1, z = value,fill = value))
  v2 <- fig + geom_tile() +  geom_contour(colour = "white")  + 
    scale_fill_gradientn(colours = heat.colors(500), name='prob') 
  v2 + xlab("n") + ylab("alpha")


# figure for the entire rank consistency
  set.seed(1)
  x = matrix(rnorm(3e+6),,3)
  n = 125
  prob = c()
  for (iter in 1:length(gamma.vec))
  {
    pop_min = pop_list[[iter]]
    m = pop_min$est[-4]
    v = pop_min$cov/n
    fit = svd(v)
    #A matrix
    #A = fit$u%*%diag(sqrt(fit$d))
    #t(A) matrix
    At = diag(sqrt(fit$d))%*%t(fit$u)
    Ax = x%*%At
    Ax = Ax + rep(m, each = nrow(Ax))
    idx =(Ax[,1] > Ax[,2]) &  (Ax[,2] > Ax[,3]) & (Ax[,3]>0)
    prob[iter] = mean(idx)
  }
  plot(prob)
  
    



plot(rmat[,10]) ; abline(h = 0.5)
which.max(rmat[,10])
head(pop_beta_mat)

# image(x = gamma.vec, y = nvec, rankcons.mat, ylab = 'n', xlab = 'gamma',
#       useRaster = TRUE)

# numerical simulations of BT
n = 6
iter = 1
iter.num = 10000
beta_mat = NULL

for (iter in 1:iter.num)
{
  Gmat_hat = Qmat = matrix(0,p,p)
  for (i in 1:(p-1) )
  {
    for (j in (i+1):p)
    {
      un = n
      if ( (i == 2 & j == 3) ) un = alpha*n
      if ( (i == 1 & j == 4) ) un = alpha*n
      if ((i==3) & (j==4)) un = alpha*n
      v = mean( ( (rt(un,df) + lambda.vec[i]) - 
                    (rt(un,df) + lambda.vec[j])  ) >0 )                
      Gmat_hat[i,j] = v 
      Gmat_hat[j,i] = 1-v 
      Qmat[i,j] = Qmat[j,i] = un
    }
  }
  # set an intial parameter for btFUN
  Qpmat = Qmat/sum(Qmat)*2
  Qmat_fit$Qpmat = Qpmat
  Qmat_fit$Gmat_hat = Gmat_hat
  Qmat_fit$x = x
  Qmat_fit$y = y
  # fit the BT
  beta_hat = btFun(Qmat_fit)  
  beta_mat = rbind(beta_mat, beta_hat)
}

sum(Qmat)/2
mean(beta_mat[,1] - beta_mat[,2]>0)
# 
# the ith row -> the result of gamma.vec[i]
# gamma.v = 0 denotes BT and gamma.v = 1 denotes gBT.
#rmat[1,]
# the jth col -> the result of nvec[j]
#rmat[,1]
j = 34
rmat[1,j]


abline(h = 0)

# numerical simulations of gBT
jj = 57
gamma.v = gamma.vec[jj]
n = 20
iter = 1
iter.num = 10000
beta_mat = NULL

for (iter in 1:iter.num)
{
  Gmat_hat = Qmat = matrix(0,p,p)
  for (i in 1:(p-1) )
  {
    for (j in (i+1):p)
    {
      un = n
      if ( (i == 2 & j == 3) ) un = alpha*n
      if ( (i == 1 & j == 4) ) un = alpha*n
      if ((i==3) & (j==4)) un = alpha*n
      v = mean( ( (rt(un,df) + lambda.vec[i]) - 
                    (rt(un,df) + lambda.vec[j])  ) >0 )                
      Gmat_hat[i,j] = v 
      Gmat_hat[j,i] = 1-v 
      Qmat[i,j] = Qmat[j,i] = un
    }
  }
  # set an intial parameter for btFUN
  Qpmat = Qmat/sum(Qmat)*2
  wmat = (1/Qpmat)^gamma.v ; diag(wmat) = 0
  Qmat_fit$Qpmat = Qpmat*wmat
  Qmat_fit$Gmat_hat = Gmat_hat
  Qmat_fit$x = x
  Qmat_fit$y = y
  # fit the BT
  beta_hat = btFun(Qmat_fit)  
  beta_mat = rbind(beta_mat, beta_hat)
}

sum(Qmat)/2
mean(beta_mat[,1] - beta_mat[,2]>0)
# 
# the ith row -> the result of gamma.vec[i]
# gamma.v = 0 denotes BT and gamma.v = 1 denotes gBT.
#rmat[1,]
# the jth col -> the result of nvec[j]
#rmat[,1]
j = which.min(abs(nvec-sum(Qmat)/2))
rmat[jj,j]


abline(h = 0)

