rm(list = ls())
setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
load("Real_BT_gBT2_cv5_all_data.rdata")
sel_idx = which(BT_est_rank <=5)
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
#set.seed(1)
vmat = NULL
for (i in 1:50)
{
  set.seed(i)
  n = nrow(rdata[,18:33])
  s_idx = sample(1:n, trunc(n*0.7))
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1
  Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx)  
  bt_est <- btFun(Qmat_fit)
  u = sort(unique(c(Qmat_fit$Qpmat)))
  gbt_est <- gbtFun(Qmat_fit, cvec = u[1])$gbt_est
  if (is.null(gbt_est)) next
  
  race_mat <- as.matrix(rdata[-s_idx,18:33])
  Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx)  
  v1 = evalFun_1(rdata, bt_est, sel_idx)
  v2 = evalFun_1(rdata, gbt_est, sel_idx)
  
  v3 = evalFun_2(rdata, bt_est, sel_idx)
  v4 = evalFun_2(rdata, gbt_est, sel_idx)
  
  v5 = evalFun_3(Qmat_fit, bt_est)
  v6 = evalFun_3(Qmat_fit, gbt_est)
  
  vmat = rbind(vmat, 
               matrix(c(v1,v2,v3,v4,v5,v6),1,6))
  
}
boxplot(vmat[,1:2])
boxplot(vmat[,3:4])
boxplot(vmat[,5:6])

