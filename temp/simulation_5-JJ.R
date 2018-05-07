## description
## investigate the variance of the proposed estimator
rm(list = ls()); gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
require('glmnet')
library(igraph)
library(ggplot2)
library(dplyr)
source('./lib/sim.R')
################################################################################
max.k = 10
p = 10
kn <- 7  ## d
rn <- 3   ## n_s
df = 1    ## degree of freedom
counter = 1
sim.iter = 200    ## 전체 simulation 과정 반복 수 
source('./lib/exe-2.R') # return the object, dmat
tn_vec = c(500,5000,50000)
cor.naive_list = list()
cor.r_list = list()
k_fold = 5
tn_i = 2

# simulation: number of obs.
  tn = tn_vec[tn_i]  ## tn 정의 (전체 rank pair의 수.)
  cat ('total n:' , tn , '\n')
  cor.naive<- rep(0,sim.iter) ## BT를 이용한 kendall's tau 저장하는 벡터 
  
  
  gbt_trueMat <- matrix(NA,sim.iter,max.k)
  gbt_truevec <- rep(NA, sim.iter)
  bt_truevec <- rep(NA, sim.iter)
  
  cor.cv <- matrix(0,sim.iter,max.k)
  cor.cv.list = vector(mode = 'list', length = sim.iter)
  ii = 1
  for  (ii in 1:sim.iter)
  {
    cat(' ',ii,'-th iteration\n')
    set.seed(ii+123) ## set seed
    Qmat = sparse_gen_fun(dmat, kn, rn, tn)
    gen_fit = gen_sim_fun(Gmat, Qmat)
    cvec<-(0:(max.k-1))/tn  ## cvec : threshold cval
    
    # cross validation
    fit = cv.gbt_fun(gen_fit, cvec, k_fold, lambda.vec)
    cv_mean_vec = colMeans(fit, na.rm = TRUE)
    cor.cv.list[[ii]] = fit
    cor.cv[ii,] =  cv_mean_vec
    
    # cor.r
    for (k in 1:length(cvec))
    {
      cval = cvec[k]
      gbt_fit = gbt_fun(gen_fit, cval, lambda.vec)
      gbt_trueMat[ii,k] = gbt_fit$cor
    }
    
    # cv selection
    min_idx <-which.max(cv_mean_vec)
    cval <-  cvec[min_idx]

    gbt_fit = gbt_fun(gen_fit, cval, lambda.vec)
    bt_fit = bt_fun(gen_fit, lambda.vec)
    
    gbt_truevec[ii] = gbt_fit$cor
    bt_truevec[ii] = bt_fit$cor
  } 
boxplot(gbt_truevec,bt_truevec)      
save.image("sim_result_5-2.rdata")      

colMeans(cor.r_list1[[1]] , na.rm=T)
colMeans(cor.r_list1[[2]] , na.rm=T)
colMeans(cor.r_list1[[3]] , na.rm=T)

colMeans(cor.r_list3[[1]] , na.rm=T)
colMeans(cor.r_list3[[2]] , na.rm=T)
colMeans(cor.r_list3[[3]] , na.rm=T)

colMeans(cor.r_list5[[1]] , na.rm=T)
colMeans(cor.r_list5[[2]] , na.rm=T)
colMeans(cor.r_list5[[3]] , na.rm=T)

setwd("E:\\rank consistency\\simulation\\Rdata\\")
save.image('.\\Simulation2_dh_170613.rdata')    

## plot 그리기 
setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim5_dh_ver2_1.jpg')
plot(0:9 , colMeans(cor.r_list1[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
lines(0:9 , colMeans(cor.r_list1[[2]] , na.rm=T) , type='b' , lty=2)
lines(0:9 , colMeans(cor.r_list1[[3]] , na.rm=T) , type='b' , lty=3)
abline(v=1 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver2_3.jpg')
plot(0:9 , colMeans(cor.r_list3[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
lines(0:9 , colMeans(cor.r_list3[[2]] , na.rm=T) , type='b' , lty=2)
lines(0:9 , colMeans(cor.r_list3[[3]] , na.rm=T) , type='b' , lty=3)
abline(v=3 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver2_5.jpg')
plot(0:9 , colMeans(cor.r_list5[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
lines(0:9 , colMeans(cor.r_list5[[2]] , na.rm=T) , type='b' , lty=2)
lines(0:9 , colMeans(cor.r_list5[[3]] , na.rm=T) , type='b' , lty=3)
abline(v=5 , col='red')
dev.off()


################################################################################
## 하나씩 따로 그리기
setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim5_dh_ver3_1_500.jpg')
plot(0:9 , colMeans(cor.r_list1[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=1 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_1_5000.jpg')
plot(0:9 , colMeans(cor.r_list1[[2]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=1 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_1_50000.jpg')
plot(0:9 , colMeans(cor.r_list1[[3]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=1 , col='red')
dev.off()


jpeg('.\\sim5_dh_ver3_3_500.jpg')
plot(0:9 , colMeans(cor.r_list3[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation' , 
     cex.lab=1.8, cex.axis=1.5, cex.sub=2.3 , lwd=1.5)
abline(v=3 , col='red' , lwd=1.5)
dev.off()

jpeg('.\\sim5_dh_ver3_3_5000.jpg')
plot(0:9 , colMeans(cor.r_list3[[2]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation' , 
     cex.lab=1.8, cex.axis=1.5, cex.sub=2.3 , lwd=1.5)
abline(v=3 , col='red' , lwd=1.5)
dev.off()

jpeg('.\\sim5_dh_ver3_3_50000.jpg')
plot(0:9 , colMeans(cor.r_list3[[3]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation' , 
     cex.lab=1.8, cex.axis=1.5, cex.sub=2.3 , lwd=1.5)
abline(v=3 , col='red' , lwd=1.5)
dev.off()


jpeg('.\\sim5_dh_ver3_5_500.jpg')
plot(0:9 , colMeans(cor.r_list5[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=5 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_5_5000.jpg')
plot(0:9 , colMeans(cor.r_list5[[2]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=5 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_5_50000.jpg')
plot(0:9 , colMeans(cor.r_list5[[3]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=5 , col='red')
dev.off()

