# the step of the experiment
# 1. select the top-13 ranked cars and fit the ranking 
# 2. select the all cars and fit the ranking 
# 3. select the all cars and fit the ranking without the result of the top-13 cars
# intransitivity


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
n = nrow(rdata)

bt_est.list = gbt_est.list = sr1_est.list =
  sr_est.list = gbt_est.list2 =
  vector(mode='list', length = sim.num)

bt_result.list = gbt_result.list = sr1_result.list =
  sr_result.list = gbt_result.list2 =
  vector(mode='list', length = sim.num)
# training code
race_mat <- as.matrix(rdata[,18:33])
num_vec <- rdata$V1
Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                   p=43, sel_idx)  
# 1:13
bt_est <- rank(btFun(Qmat_fit))
bt_result <- make_result(bt_est)
sr_est <- srFun(Qmat_fit)
sr_result  <- make_result(sr_est)
sr1.result = sr1_fun(Qmat_fit)
sr1_est = gbtFun_recov(sr1.result, Qmat_fit, method='count',
                       allowties = F)

Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                   p=43, sel_idx)  
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
gbt_fit_result = gbt_fit$sc_list
gbt_est = gbtFun_recov(gbt_fit_result, Qmat_fit, 
                       method = 'count', allowties = F)





evalFun_3_pair(bt_result, Qmat_fit)
evalFun_3_pair(gbt_fit_result, Qmat_fit)
evalFun_3_pair(sr1_fun(Qmat_fit), Qmat_fit)

evalFun_3(Qmat_fit, bt_est)
evalFun_3(Qmat_fit, gbt_est)
evalFun_3(Qmat_fit, sr1_est)
evalFun_3(Qmat_fit, sr_est)

## off_set 
Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                   p=43, sel_idx = setdiff(sel_idx,40), 
                   off_set = T)  
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
gbt_est_off = gbtFun_recov(gbt_fit$sc_list, Qmat_fit, 
                       method = 'count', allowties = F)
gbt_est_off = rank(gbt_est_off[sel_idx])
evalFun_3(Qmat_fit, gbt_est_off)

bt_est_off <- btFun(Qmat_fit)
bt_est_off = rank(bt_est_off[sel_idx])
evalFun_3(Qmat_fit, bt_est_off)


sr_est_off <- srFun(Qmat_fit)
sr_est_off = rank(sr_est_off[sel_idx])
evalFun_3(Qmat_fit, sr_est_off)

