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

## off_set 
Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                   p=43, sel_idx = setdiff(sel_idx,40), 
                   off_set = T)  
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
gbt_est_off = gbtFun_recov(gbt_fit_result, Qmat_fit, 
                       method = 'count', allowties = F)
gbt_est_off = rank(gbt_est_off[sel_idx])


evalFun_3(Qmat_fit, gbt_est_off)


#1:43
Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                   p=43, sel_idx = 1:43)  
bt_est43 <- rank(btFun(Qmat_fit))
bt_est43  = rank(bt_est43[sel_idx])

#bt_result <- make_result(bt_est)
sr1.result = sr1_fun(Qmat_fit)
sr1_est43 = gbtFun_recov(sr1.result, Qmat_fit, method='count',
                       allowties = F)
sr1_est43 = rank(sr1_est43[sel_idx])

Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                   p=43, sel_idx = 1:43)  
gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
gbt_fit_result = gbt_fit$sc_list
gbt_est43 = gbtFun_recov(gbt_fit_result, Qmat_fit, 
                       method = 'count', allowties = F)
gbt_est43 = rank(gbt_est43[sel_idx])

cbind(bt_est, bt_est43)
cbind(sr1_est, sr1_est43)
cbind(gbt_est, gbt_est43)


evalFun_3_pair(bt_result, Qmat_fit, sel_idx)
evalFun_3_pair(gbt_fit_result, Qmat_fit, sel_idx)
evalFun_3_pair(sr1_fun(Qmat_fit), Qmat_fit, sel_idx )



a1 = evalFun_4_pair(bt_result, Qmat_fit, sel_idx)
a2 = evalFun_4_pair(gbt_fit_result, Qmat_fit, sel_idx)  
a3 = evalFun_4_pair(sr1_fun(Qmat_fit), Qmat_fit, sel_idx)
cbind(a1$v1,a2$v1,a3$v1)
#evalFun_3_pair(sr2_fun(Qmat_fit), Qmat_fit)



# evalFun_4_pair(bt_result, Qmat_fit)
# evalFun_4_pair(gbt_fit_result, Qmat_fit)  
sr1_result = sr1_fun(Qmat_fit)
# evalFun_4_pair(sr1_fun(Qmat_fit), Qmat_fit)
# 


b1 = bt_result[which(bt_result[,2]==43),]
b2 = gbt_fit_result[which(gbt_fit_result[,2]==43),]
b3 = sr_result[which(sr_result[,2]==43),]

### check the result of car43
cbind(b1[,1:3],b2[,3],b3[,3])

# car 43 vs 4, 7, 8, 12, 13, 19, 20, 38

# car 43 vs car 4
