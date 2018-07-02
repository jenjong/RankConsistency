rm(list = ls())
setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
load("Real_BT_gBT2_cv5_all_data.rdata")
gBT2_est_rank
BT_est_rank
# 
which(BT_est_rank==9)
gBT2_est_rank[40]
which(BT_est_rank==4)
gBT2_est_rank[12]
sel_idx = which(BT_est_rank <=13)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')

rdata<-read.csv('racing_data.csv', header=F)
rdata <- rbind(rdata, rdata)
max_k = 4
cvec_r <- seq(0, max_k, by = 2)
file_idx = 1
cat("iteration::", seed_v, '\n')
seed_v_i = (file_idx -1)*inner_iter + seed_v
set.seed(seed_v_i)
sc_list = vector(mode ='list', length = max_k)
sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*1)))  
race_mat <- as.matrix(rdata[sample_idx,18:33])   ## train set의 각 게임당 선택 차종 
num_vec<- rdata$V1[sample_idx]  ## 각 게임마다 참여한 유저 수 
Qmat_fit <-QmatFunc(race_mat, num_vec)  
Qmat = Qmat_fit$Qmat
  
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
  