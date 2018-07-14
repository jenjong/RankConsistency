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
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, ctype = 'balance')
result = gbt_fit$sc_list
gbt_est = gbtFun_recov(result, Qmat_fit, method = 'count', allowties = F)

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

sr1_est = gbtFun_recov(result, Qmat_fit, method='count', allowties = F)
evalFun_3(Qmat_fit, sr1_est)

a1 = evalFun_4(Qmat_fit, gbt_est)
a2 = evalFun_4(Qmat_fit, sr1_est)
a3 = cbind(a1$v1, a2$v1, a1$v2)
colSums(a3)/98
colMeans(a3[,1:2]/a3[,3])

result1 = gbt_fit$sc_list
a1 = evalFun_4_pair(result1, Qmat_fit)
result2 = sr1_fun(Qmat_fit)
a2 = evalFun_4_pair(result2, Qmat_fit)
a3 = cbind(a1$v1, a2$v1, a1$v2)

# 
gbt_bt = evalFun_5(Qmat_fit, gbt_est)$v1
sr_bt = evalFun_5(Qmat_fit, sr1_est)$v1
#
gbt = evalFun_5_pair(result1, Qmat_fit)$v1
sr = evalFun_5_pair(result2, Qmat_fit)$v1
Qmat = evalFun_5_pair(result2, Qmat_fit)$v2

Nmat = Qmat_fit$Qmat
save(list = c("gbt_bt", "sr_bt", "gbt", "sr", "Qmat", "Nmat"),
     file = "real_sim.rdata")

a2 = evalFun_4(Qmat_fit, sr1_est)
a3 = cbind(a1$v1, a2$v1, a1$v2)




result = gbt_fit$sc_list
gbt_est = gbtFun_recov(result, Qmat_fit, method = 'count')
evalFun_3(Qmat_fit, gbt_est)


