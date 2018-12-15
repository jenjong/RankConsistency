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
avec = seq(0,0.4, by = 0.05)
for (i in 1:sim.num)
{
  cat(i,'\n')
  set.seed(i)
  s_idx = sample(1:n, trunc(n*0.7))
  # training code
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1[s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  race_mat <- as.matrix(rdata[-s_idx,18:33])
  num_vec <- rdata$V1[-s_idx]
  Qmat_fit_ts <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  
  
  #bt_est <- btFun(Qmat_fit)
  # a = 0, gbt is bt, a  = 1, gbt is full-weight bt
  v = c()
  for (k in 1:length(avec))
  {
    gbt_est <- gbtFun_ver01(Qmat_fit, a = avec[k])
    v[k] = evalFun_3(Qmat_fit_ts, gbt_est)
  }
  vmat = rbind(vmat, v)
  #bt_result <- make_result(bt_est)
  #gbt_result <- make_result(gbt_est)
  #sr_est <- srFun(Qmat_fit)
  #sr_result  <- make_result(sr_est)
  
  #v1 = evalFun_3(Qmat_fit_ts, bt_est)
  #vmat = rbind(vmat, c(v1,v2,v3))
}
colMeans(vmat)

