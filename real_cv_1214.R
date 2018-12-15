rm(list = ls())
gc()
# training code
# set path
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/GitHub/RankConsistency")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
}
# load car segmentation
library(igraph)
library(MASS)
library(glmnet)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
sim.num = 100

rdata<-read.csv('racing_data.csv', header=F)
sel_idx = 1:43
n = nrow(rdata)

bt_est.list = gbt_est.list = sr1_est.list =
  sr_est.list =   vector(mode='list', length = sim.num)

bt_result.list = gbt_result.list = sr1_result.list =
  sr_result.list  =
  vector(mode='list', length = sim.num)
i = 1
vmat = NULL
gbt_est_mat = NULL
avec = seq(0,0.3, by = 0.05)
sel.a = c()
for (i in 101:200)
{
  cat(i,'\n')
  set.seed(i)
  s_idx = sample(1:n, trunc(n*0.7))
  cv_idx = sample(1:5, length(s_idx), replace = T)
  # training code
  j = 1
  v = rep(0, length(avec))
  for (j in 1:5)
  {
    cv_tr_idx = s_idx[cv_idx!=j]
    cv_ts_idx = s_idx[cv_idx==j]
    race_mat <- as.matrix(rdata[cv_tr_idx,18:33])
    num_vec <- rdata$V1[cv_tr_idx]
    Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                       p=43, sel_idx)  
    race_mat <- as.matrix(rdata[cv_ts_idx,18:33])
    num_vec <- rdata$V1[cv_ts_idx]
    Qmat_fit_ts <-QmatFun(race_mat, num_vec, cut_var = 0,
                          p=43, sel_idx)
    for (k in 1:length(avec))
    {
      gbt_est <- gbtFun_ver01(Qmat_fit, a = avec[k])
      v[k] = v[k] + evalFun_3(Qmat_fit_ts, gbt_est)
    }
  }
  
  k <- which.max(v)
  sel.a[i] = k
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1[s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  race_mat <- as.matrix(rdata[-s_idx,18:33])
  num_vec <- rdata$V1[-s_idx]
  Qmat_fit_ts <-QmatFun(race_mat, num_vec, cut_var = 0,
                        p=43, sel_idx)
  bt_est <- gbtFun_ver01(Qmat_fit, a = 0)
  gbt_est <- gbtFun_ver01(Qmat_fit, a = avec[k])
  v1 = evalFun_3(Qmat_fit_ts, bt_est)
  v2 = evalFun_3(Qmat_fit_ts, gbt_est)
  vmat = rbind(vmat, c(v1,v2))
  cat("perf::", colMeans(vmat), "\n")
  cat("selection::", table(sel.a), "\n")
}
  
colMeans(vmat)

