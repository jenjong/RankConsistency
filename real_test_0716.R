rm(list = ls())
gc()
# training code
# set path
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/Github/RankConsistency")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
}
# load car segmentation
load("Real_BT_gBT2_cv5_all_data.rdata")
i_1 = 1
i_2 = 13
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)

# library 
library(MASS)
library(igraph)
library(glmnet)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
sim.num = 50

rdata<-read.csv('racing_data.csv', header=F)
rdata = rbind(rdata,rdata)
n = nrow(rdata)

load(paste('./result/real_traninig', i_1, i_2, sep='_'))
# test procedure
vmat1 = vmat2 = vmat3 = NULL
i = 1
for (i in 1:sim.num)
{
  set.seed(i)
  s_idx = sample(1:n, trunc(n*0.7))
  
  race_mat <- as.matrix(rdata[-s_idx,18:33])
  num_vec <- rdata$V1[-s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  
  gbt_est = gbt_est.list[[i]]
  if (is.null(gbt_est)) next
  bt_est = bt_est.list[[i]]
  sr_est = sr_est.list[[i]]
  sr1_est = sr1_est.list[[i]]
  
  v1 = evalFun_1(rdata[-s_idx,], bt_est, sel_idx)
  v2 = evalFun_1(rdata[-s_idx,], gbt_est, sel_idx)
  v3 = evalFun_1(rdata[-s_idx,], sr1_est, sel_idx)
  v4 = evalFun_1(rdata[-s_idx,], sr_est, sel_idx)
  vmat1 = rbind(vmat1, c(v1,v2,v3,v4))
  
  v1 = evalFun_2(rdata[-s_idx,], bt_est, sel_idx)
  v2 = evalFun_2(rdata[-s_idx,], gbt_est, sel_idx)
  v3 = evalFun_2(rdata[-s_idx,], sr1_est, sel_idx)
  v4 = evalFun_2(rdata[-s_idx,], sr_est, sel_idx)
  vmat2 = rbind(vmat2, c(v1,v2,v3,v4))

  v1 = evalFun_3(Qmat_fit, bt_est)
  v2 = evalFun_3(Qmat_fit, gbt_est)
  v3 = evalFun_3(Qmat_fit, sr1_est)
  v4 = evalFun_3(Qmat_fit, sr_est)
  vmat3 = rbind(vmat3, c(v1,v2,v3,v4))
}

boxplot(vmat1[,1:4], names= c("a",'b','c', 'd'))
boxplot(vmat2[,1:4], names= c("a",'b','c', 'd'))
boxplot(vmat3[,1:4], names= c("a",'b','c', 'd'))


colMeans(vmat1[,1:3])
colMeans(vmat2[,1:3])
colMeans(vmat3[,1:3])
