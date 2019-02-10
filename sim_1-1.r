rm (list =ls()); gc()
library(e1071)
library(quadprog)
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
source("./lib/lib_rank.r")
#source('D:\\JJJ\\ranking problem\\Jeon-2014-03-05 order consistency\\experiment\\library\\lib_rank.r')

    n=200 ; df = 1
    lambda.vec = c(1.2,1.1,0.6,0)*4
    p = length(lambda.vec)
    alpha.vec = 1:20
    num.iter = 100
    beta.mat1 = beta.mat2  = NULL
    rank.mat1 = rank.mat2 = NULL 
    set.seed(1)
    for (k in 1:length(alpha.vec))
    {
    alpha = alpha.vec[k]
    tmp.mat1 = tmp.mat2 = NULL
        for ( iter in 1:100)
        {
        X = NULL
        pij = c()
        true.prob = c()
        idx = 1

              for (i in 1:(p-1) )
              {
                    for (j in (i+1):p)
                    {
                    un = n
                    if ( (i == 2 & j == 3) ) un = n*alpha
                    if ( (i == 1 & j == 4) ) un = n*alpha
                    if ((i==3) & (j==4)) un=n*alpha
                    tmp.vec = rep(0,p)
                    tmp.vec[i] = 1 ; tmp.vec[j] = -1
                    X = rbind(X, tmp.vec)
                    true.prob[idx] = mean( ( (rt(un,df) + lambda.vec[i]) - (rt(un,df) + lambda.vec[j])  ) >0 )
                    #          true.prob[idx]=pt(lambda.vec[i]-lambda.vec[j],1)
                    pij[idx] = un
                    idx = idx + 1
                    }
              }
        # orginal Bradley-Terry model
        X = X[,-p] ; pij = pij/ sum(pij)
        beta.vec = IWLS.ridge(X,true.prob,pij,accuracy=10e-10,maxitration=200)
        beta.vec = c(beta.vec,0)
        tmp.mat1 = rbind(tmp.mat1,beta.vec)

        beta.vec = pIWLS(X,true.prob,pij,accuracy=10e-10,maxitration=200)
        tmp.mat2 = rbind(tmp.mat2,beta.vec)

#        rownames(X) = NULL
#        m = svm(mX,my, kernel='linear', scale = F, cost = 0.5)
#        w <- t(m$coefs) %*% X[m$index,]
#        w = c(w,0)
#        w
#        cat(beta.vec,'\n')
        }

    beta.mat1 = rbind(beta.mat1, colMeans(tmp.mat1))
    beta.mat2 = rbind(beta.mat2, colMeans(tmp.mat2))
    rank.mat1 = rbind(rank.mat1, colMeans(rank.fun(tmp.mat1,p)))
    rank.mat2 = rbind(rank.mat2, colMeans(rank.fun(tmp.mat2,p)))
    }

  # setwd('D:\\JJJ\\ranking problem\\Jeon-2014-03-05 order consistency\\rank-consistenmcy_ver04-JJ')
  # jpeg(filename = "Rplot1.jpg")
  # plot(alpha.vec, beta.mat1[,1], type = 'b', lty = 1, ylim = c(0,3),
  #      ylab = expression(beta), xlab = expression(alpha), main = 'mean of estimated coefficients',
  #      pch = 1)
  # lines(alpha.vec, beta.mat1[,2], lty = 2,type = 'b', pch = 2)
  # lines(alpha.vec, beta.mat1[,3], lty = 3, type = 'b',pch = 3)
  # lines(alpha.vec, beta.mat1[,4], lty = 4, type = 'b',pch = 4)
  # dev.off()
  # jpeg(filename = "Rplot2.jpg")
  # plot(alpha.vec, beta.mat2[,1], type = 'b', lty = 1, ylim = c(0,3),
  #      ylab = expression(beta), xlab = expression(alpha), main = 'mean of estimated coefficients',
  #      pch = 1)
  # lines(alpha.vec, beta.mat2[,2], lty = 2,type = 'b', pch = 2)
  # lines(alpha.vec, beta.mat2[,3], lty = 3, type = 'b',pch = 3)
  # lines(alpha.vec, beta.mat2[,4], lty = 4, type = 'b',pch = 4)
  # dev.off()
    
  png(filename = "Rplot3.png")
  plot(alpha.vec, rank.mat1[,1], type = 'b', lty = 1, ylim = c(0,4.2),
       ylab = "averages of estimated ranks", xlab = expression(gamma), 
       main = 'BT',
       pch = 20, cex.lab = 1.5, cex.axis = 1.3)
  lines(alpha.vec, rank.mat1[,2], lty = 1,type = 'b', pch = 2)
  lines(alpha.vec, rank.mat1[,3], lty = 1, type = 'b',pch = 3)
  lines(alpha.vec, rank.mat1[,4], lty = 1, type = 'b',pch = 4)
  legend('bottomright', legend = c(expression(paste(lambda[1], "     ")),
                                   expression(lambda[2]),
                                   expression(lambda[3]),
                                   expression(lambda[4])),
         pch = c(20, 2,3,4), lty = c(1,1,1,1) )
  
  
  dev.off()
  
  png(filename = "Rplot4.png")
  plot(alpha.vec, rank.mat2[,1], type = 'b', lty = 1, ylim = c(0,4.2),
       ylab = "averages of estimated ranks", xlab = expression(gamma), 
       main = 'gBT',
       pch = 20, cex.lab = 1.5, cex.axis = 1.3)
  lines(alpha.vec, rank.mat2[,2], lty = 1,type = 'b', pch = 2)
  lines(alpha.vec, rank.mat2[,3], lty = 1, type = 'b',pch = 3)
  lines(alpha.vec, rank.mat2[,4], lty = 1, type = 'b',pch = 4)
  legend('bottomright', legend = c(expression(paste(lambda[1], "     ")),
                                   expression(lambda[2]),
                                   expression(lambda[3]),
                                   expression(lambda[4])),
         pch = c(20, 2,3,4), lty = c(1,1,1,1) )
  
  dev.off()
