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
sim.num = 1

rdata<-read.csv('racing_data.csv', header=F)
#rdata = rbind(rdata,rdata)
n = nrow(rdata)

bt_est.list = gbt_est.list = sr1_est.list =
  sr_est.list = 
  vector(mode='list', length = sim.num)
i = 1

for (i in 1:sim.num)
{
  set.seed(i)
  s_idx = sample(1:n, trunc(n*0.7))
  
  # training code
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1[s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  
  bt_est <- btFun(Qmat_fit)
  sr_est <- srFun(Qmat_fit)
  gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'balance')
  if (is.null(gbt_fit$gbt_est)) next
  gbt_fit.result = gbt_fit$sc_list
  gbt_est = gbtFun_recov(gbt_fit.result, Qmat_fit, 
                         method = 'count', allowties = F)
  sr1.result = sr1_fun(Qmat_fit)
  sr1_est = gbtFun_recov(sr1.result, Qmat_fit, method='count',
                         allowties = F)
  sr_est.list[[i]] = sr_est
  bt_est.list[[i]] = bt_est
  gbt_est.list[[i]] = gbt_est
  sr1_est.list[[i]] = sr1_est
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
              'gbt_est.list', 'sr1_est.list'))
