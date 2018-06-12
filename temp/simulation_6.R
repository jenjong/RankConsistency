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
tn_vec = c(500,5000,100000)
cor.naive_list = list()
cor.r_list = list()
k_fold = 5
tn_i = 3
cor.cv <- matrix(0,sim.iter,3)
for (tn_i in 1:3)
{
  # simulation: number of obs.
  tn = tn_vec[tn_i]  ## tn 정의 (전체 rank pair의 수.)
  cat ('total n:' , tn , '\n')
  cor.naive<- rep(0,sim.iter) ## BT를 이용한 kendall's tau 저장하는 벡터 
  
  
  gbt_trueMat <- matrix(NA,sim.iter,max.k)
  gbt_truevec <- rep(NA, sim.iter)
  bt_truevec <- rep(NA, sim.iter)
  

  ii = 1
  for  (ii in 1:sim.iter)
  {
    cat(' ',ii,'-th iteration\n')
    set.seed(ii+123) ## set seed
    Qmat = sparse_gen_fun(dmat, kn, rn, tn)
    gen_fit = gen_sim_fun(Gmat, Qmat)
    Gmat.hat <- gen_fit$G
    Qmat <- gen_fit$Q
    p = ncol(Qmat)
    Gmat.hat <- Gmat.hat/Qmat
    Gmat.hat[!is.finite(Gmat.hat)] = 0.5
    diag(Gmat.hat) = 0
    cor.cv[ii, tn_i] = cor(p:1, apply(Gmat.hat, 1, sum), method = 'kendall')
  }
}

  boxplot(cor.cv[,3])
  
  mean(cor.cv[,1])
  
  ###### real
  