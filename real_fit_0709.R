rm(list = ls())
setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
load("Real_BT_gBT2_cv5_all_data.rdata")
gBT2_est_rank
BT_est_rank
#which(BT_est_rank==9)
#gBT2_est_rank[40]
#which(BT_est_rank==4)
#gBT2_est_rank[12]
sel_idx = which(BT_est_rank <=13)
#sel_idx = sel_idx[!sel_idx==8]
#sel_idx = sel_idx[!sel_idx==20]
#sel_idx = sel_idx[!sel_idx==34]
#sel_idx = sel_idx[!sel_idx==38]
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')
rdata<-read.csv('racing_data.csv', header=F)
# data preprocessing
race_mat <- as.matrix(rdata[,18:33])
num_vec<- rdata$V1
Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx)  


# estimation
bt_est <- btFun(Qmat_fit)
u = sort(unique(c(Qmat_fit$Qpmat)))
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, ctype = 'boost')
gbt_est <- gbt_fit$gbt_est

# evaluation
evalFun_1(rdata, bt_est, sel_idx)
evalFun_1(rdata, gbt_est, sel_idx)

evalFun_2(rdata, bt_est, sel_idx)
evalFun_2(rdata, gbt_est, sel_idx)


evalFun_3_pair(gbt_fit$sc_list, Qmat_fit)
evalFun_3_pair(sr1_fun(Qmat_fit), Qmat_fit)
evalFun_3_pair(sr2_fun(Qmat_fit), Qmat_fit)


evalFun_3(Qmat_fit, bt_est)
evalFun_3(Qmat_fit, gbt_est)
result = sr1_fun(Qmat_fit)

sr1_est = gbtFun_recov(result, Qmat_fit, method='gaussian')
evalFun_3(Qmat_fit, sr1_est)

result = gbt_fit$sc_list
gbt_est = gbtFun_recov(result, Qmat_fit, method = 'gaussian')
evalFun_3(Qmat_fit, gbt_est)

# 
set.seed(1)
n = nrow(rdata[,18:33])
s_idx = sample(1:n, trunc(n*0.7))
race_mat <- as.matrix(rdata[s_idx,18:33])
num_vec <- rdata$V1
Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx)  
bt_est <- btFun(Qmat_fit)
u = sort(unique(c(Qmat_fit$Qpmat)))
gbt_est <- gbtFun(Qmat_fit, cvec = u[2])$gbt_est


##


plot(rank(bt_est), rank(gbt_est))


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