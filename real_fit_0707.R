rm(list = ls())
setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
load("Real_BT_gBT2_cv5_all_data.rdata")
sel_idx = which(BT_est_rank <=13)
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
library('glmnet')
rdata<-read.csv('racing_data.csv', header=F)
vmat = NULL
i=1
for (i in 1:50)
{
  set.seed(i)
  n = nrow(rdata[,18:33])
  s_idx = sample(1:n, trunc(n*0.7))
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1[s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                     p=43, sel_idx)  
  bt_est <- btFun(Qmat_fit)
  gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
  if (is.null(gbt_fit$gbt_est)) next
  
  gbt_fit.result = gbt_fit$sc_list
  gbt_est = gbtFun_recov(gbt_fit.result, Qmat_fit, 
                         method = 'count', allowties = F)
  sr1.result = sr1_fun(Qmat_fit)
  sr1_est = gbtFun_recov(sr1.result, Qmat_fit, method='count',
                         allowties = F)
  race_mat <- as.matrix(rdata[-s_idx,18:33])
  num_vec <- rdata$V1[-s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  # v1 = evalFun_1(rdata, bt_est, sel_idx)
  # v2 = evalFun_1(rdata, gbt_est, sel_idx)
  # v3 = evalFun_1(rdata, sr1_est, sel_idx)
  # vmat = rbind(vmat, c(v1,v2,v3))

  v1 = evalFun_3(Qmat_fit, bt_est)
  v2 = evalFun_3(Qmat_fit, gbt_est)
  v3 = evalFun_3(Qmat_fit, sr1_est)
  vmat = rbind(vmat, c(v1,v2,v3))
}

evalFun_3(Qmat_fit, bt_est)
evalFun_3(Qmat_fit, gbt_est)
result = sr1_fun(Qmat_fit)
sr1_est = gbtFun_recov(result, Qmat_fit, method='count', allowties = F)
evalFun_3(Qmat_fit, sr1_est)

boxplot(vmat[,1:3])
colMeans(vmat[-10,])
boxplot(vmat[-10,1:3])

s_idx = sample(1:n, trunc(n*1))
race_mat <- as.matrix(rdata[s_idx,18:33])
num_vec <- rdata$V1[s_idx]
Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx)  
bt_est <- btFun(Qmat_fit)
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'none')
gbt_est <- gbt_fit$gbt_est
6- rank(bt_est)
6- rank(gbt_est)
