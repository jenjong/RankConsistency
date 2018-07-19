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
i_2 = 43
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
#rdata = rbind(rdata,rdata)
n = nrow(rdata)

bt_est.list = gbt_est.list = sr1_est.list =
  sr_est.list = gbt_est.list2 =
  vector(mode='list', length = sim.num)

bt_result.list = gbt_result.list = sr1_result.list =
  sr_result.list = gbt_result.list2 =
  vector(mode='list', length = sim.num)
i = 1

for (i in 1:sim.num)
{
  cat(i,'\n')
  set.seed(i)
  s_idx = sample(1:n, trunc(n))
  
  # training code
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1[s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  
  bt_est <- btFun(Qmat_fit)
  bt_result <- make_result(bt_est)
  sr_est <- srFun(Qmat_fit)
  sr_result  <- make_result(sr_est)
  sr1.result = sr1_fun(Qmat_fit)
  sr1_est = gbtFun_recov(sr1.result, Qmat_fit, method='count',
                         allowties = F)

  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                     p=43, sel_idx)  
  gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
  if (is.null(gbt_fit$gbt_est)) next
  gbt_fit.result = gbt_fit$sc_list
  gbt_est = gbtFun_recov(gbt_fit.result, Qmat_fit, 
                         method = 'count', allowties = F)

  sr_result.list[[i]] = sr_result
  bt_result.list[[i]] = bt_result
  gbt_result.list[[i]] = gbt_fit.result
  sr1_result.list[[i]] = sr1.result
  
  
  sr_est.list[[i]] = sr_est
  bt_est.list[[i]] = bt_est
  gbt_est.list[[i]] = gbt_est
  sr1_est.list[[i]] = sr1_est

  gbt_est_bin = gbtFun_recov(gbt_fit.result, Qmat_fit, 
                         method = 'binomial', allowties = F)
  gbt_est.list2[[i]] = gbt_est_bin
  #if (any(gbt_est_bin != gbt_est)) break
}

if (Sys.info()[1] == "Linux")
{
  restorePath = '/home/jeon/Dropbox/GitHub/RankConsistency'
} else {
  restorePath = 'C:/Users/Jeon/Dropbox/GitHub/RankConsistency'
}

save(file = paste0(restorePath,
                  '/result/real_traninig_', i_1,"_", i_2), 
     list = c('sr_est.list',   'bt_est.list', 
              'gbt_est.list', 'sr1_est.list',
              'gbt_est.list2'))
