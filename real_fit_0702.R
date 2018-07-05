rm(list = ls())
setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
load("Real_BT_gBT2_cv5_all_data.rdata")
gBT2_est_rank
BT_est_rank
# 
#which(BT_est_rank==9)
#gBT2_est_rank[40]
#which(BT_est_rank==4)
#gBT2_est_rank[12]
sel_idx = which(BT_est_rank <=13)
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')
rdata<-read.csv('racing_data.csv', header=F)
max_k = 0
cvec_r <- 0
sc_list = vector(mode ='list', length = max_k)
# data preprocessing
race_mat <- as.matrix(rdata[,18:33])
num_vec<- rdata$V1
Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx)  
bt_est <- btFun(Qmat_fit)
gbt_est <- gbtFun(Qmat_fit)$gbt_estmat
evalFun_1(rdata, bt_est, sel_idx)
evalFun_1(rdata, gbt_est, sel_idx)




  sum(Qmat[12,-sel_idx])
  sum(Qmat[12,sel_idx])
  sum(Qmat[40,-sel_idx])
  sum(Qmat[40,])
  Gmat_hat <- Qmat_fit$Gmat_hat
  Wmat <- Qmat_fit$Wmat
  Gmat_hat[40,sel_idx]
  Gmat_hat[12,sel_idx]
  cbind(Gmat_hat[40,sel_idx], Gmat_hat[12,sel_idx])
  cbind(Wmat[40,sel_idx], Wmat[12,sel_idx])
  cbind(Wmat[40,sel_idx], Wmat[12,sel_idx])
  cbind(Qmat[40,sel_idx], Qmat[12,sel_idx])
  cbind(Wmat[40,sel_idx], Wmat[12,sel_idx])/
  cbind(Qmat[40,sel_idx], Qmat[12,sel_idx])
  names(Qpmat)
  
  sum(Wmat[40,sel_idx])/sum(Qmat[40,sel_idx])
  sum(Wmat[12,sel_idx])/sum(Qmat[12,sel_idx])
  
  sum(Wmat[12,])/sum(Qmat[12,])
### measure  