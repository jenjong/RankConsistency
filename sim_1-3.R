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

# set parameters in simulation
    df = 1 ;     alpha = 10
    lambda.vec = c(1.2,1.1,0.6,0)*4
# population minimizer    
    p = length(lambda.vec)
    X = NULL
    pij = c()
    true.prob = c()
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
          X = rbind(X, tmp.vec)
          v = pt((lambda.vec[i]- lambda.vec[j])/2, df = 1)            
          true.prob[idx] = v
          Gmat[i,j] = v 
          Gmat[j,i] = 1-v 
          Qmat[i,j] = Qmat[j,i] = un
          pij[idx] = un
          idx = idx + 1
          }
    }
        # Bradley-Terry model (old version)
        # X = X[,-p] ; pij = pij/ sum(pij)
        # beta.vec = IWLS.ridge(X,true.prob,pij,accuracy=10e-10,maxitration=200)
        # beta.vec = c(beta.vec,0)
        
        # Bradley-Terry model (new version)
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
            
            # set up Qmat_fit
            Qmat_fit = list()
            Qmat_fit$Qpmat = Qpmat
            Qmat_fit$Gmat_hat = Gmat
            Qmat_fit$x = x
            Qmat_fit$y = y
            
            # fit the BT
            beta_tilde = btFun(Qmat_fit)
            beta_tilde
# BT estimator
# diagal selection:  (1 + ( 0:(p-1) )*(p+1))
# asymptotic variance
            # true probability
            true_prob = Gmat
            # variance of score function
            v_s = matrix(0,p-1,p-1)
              # diag
              for ( i in 1:(p-1))
              {
                v_s[i,i] = sum(true_prob[i,-i]*(1-true_prob[i,-i])*Qpmat[i,-i])
              }
              # off-diag
              for ( i in 1:(p-1))
              {
                for (j in 1:(p-1))
                {
                  if (i==j) next
                  v_s[i,j] = - true_prob[i,j]*(1-true_prob[i,j])*Qpmat[i,j]  
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
                H[i,i] = sum(est_prob[i,-i]*(1-est_prob[i,-i])*Qpmat[i,-i])
              }
              # off-diag
              for ( i in 1:(p-1))
              {
                for (j in 1:(p-1))
                {
                  if (i==j) next
                  H[i,j] = - est_prob[i,j]*(1-est_prob[i,j])*Qpmat[i,j]  
                }
              }
            # Asymptotic variance 
            FH = solve(H)%*%v_s%*%solve(H)
            (FH[1,1] + FH[2,2] - 2*FH[1,2])
            beta_tilde[1] - beta_tilde[2]

# gBT estimator
            # set up weight mat 
            gamma.v = 0.5
            wmat = (1/Qpmat)^gamma.v ; diag(wmat) = 0
            # set up (x,y) : complete comparisons
            # already done

            # set up Qmat_fit
            Qmat_fit = list()
            Qmat_fit$Qpmat = Qpmat*wmat
            Qmat_fit$Gmat_hat = Gmat
            Qmat_fit$x = x
            Qmat_fit$y = y
            
            # fit the gBT
            beta_tilde = btFun(Qmat_fit)        
            beta_tilde 
            
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
            (FH[1,1] + FH[2,2] - 2*FH[1,2])
            beta_tilde[1] - beta_tilde[2]
            
            
# (single) numerical simulation : compute bias for each sample size   
            # set parameters in simulation
            n = 20
            # sample generation
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
            Qpmat = Qmat/sum(Qmat)
            Qmat_fit$Qpmat = Qpmat
            Qmat_fit$Gmat_hat = Gmat_hat
            Qmat_fit$x = x
            Qmat_fit$y = y
            # fit the BT
            beta_hat = btFun(Qmat_fit)
            

                        
            
# numerical simulations : compute bias for each sample size               
            
            n = 20
            iter = 1
            iter.num = 1000
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
              Qpmat = Qmat/sum(Qmat)
              Qmat_fit$Qpmat = Qpmat
              Qmat_fit$Gmat_hat = Gmat_hat
              Qmat_fit$x = x
              Qmat_fit$y = y
              # fit the BT
              beta_hat = btFun(Qmat_fit)  
              beta_mat = rbind(beta_mat, beta_hat)
            }
            boxplot(beta_mat[,1] - beta_mat[,2])
            abline(h = 0)

            # numerical simulations : compute bias for each sample size               
            
            n = 20
            iter = 1
            iter.num = 1000
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
              Qpmat = matrix(1,p,p)
              Qmat_fit$Qpmat = Qpmat
              Qmat_fit$Gmat_hat = Gmat_hat
              Qmat_fit$x = x
              Qmat_fit$y = y
              # fit the BT
              beta_hat = btFun(Qmat_fit)  
              beta_mat = rbind(beta_mat, beta_hat)
            }
            boxplot(beta_mat[,1] - beta_mat[,2])
            abline(h = 0)            
            
            