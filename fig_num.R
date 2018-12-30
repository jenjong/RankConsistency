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
sel_idx_top = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
# library 
library(MASS)
library(igraph)
library(glmnet)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
rdata<-read.csv('racing_data.csv', header=F)
n = nrow(rdata)
s_idx = 1:n
# training code
  race_mat <- as.matrix(rdata[s_idx,18:33])
  num_vec <- rdata$V1[s_idx]
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx)  
  Gmat_hat = Qmat_fit$Gmat_hat 
  Qmat = Qmat_fit$Qmat
  Gmat_hat[Qmat == 0] = 0.5
  diag(Gmat_hat) = 0
  naive_est = apply(Gmat_hat, 1, sum)
  ar_idx = order(BT_est_rank)
  plot(naive_est[ar_idx])
  p = 43
  image(1:p,1:p,Qmat[ar_idx,ar_idx], col=rainbow(200), 
        xlab = 'j', ylab = 'k')
  abline(v = 13.5, col = 'white', lwd = 2)
  abline(h = 13.5, col = 'white', lwd = 2)
  
  image(1:p,1:p,Qmat, col=rainbow(200), 
        xlab = 'j', ylab = 'k')
  
  
  
  
  bt_est <- btFun(Qmat_fit)
  bt_result <- make_result(bt_est)
  sr_est <- srFun(Qmat_fit)
  sr_result  <- make_result(sr_est)
  sr1.result = sr1_fun(Qmat_fit)
  sr1_est = gbtFun_recov(sr1.result, Qmat_fit, method='count',
                         allowties = F)
  
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 1,
                     p=43, sel_idx)  
  gbt_fit <- gbtFun(Qmat_fit, cut_v = 0, 'none', a = 1)
  gbt_fit2 <- gbtFun(Qmat_fit, cut_v = 0, 'none', a = 0)  
  if (is.null(gbt_fit$gbt_est)) next
  gbt_fit.result = gbt_fit$sc_list
  gbt_fit.result2 = gbt_fit2$sc_list  
  gbt_est = gbtFun_recov(gbt_fit.result, Qmat_fit, 
                         method = 'count', allowties = F)
  gbt_est2 = gbtFun_recov(gbt_fit.result2, Qmat_fit, 
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
}

if (Sys.info()[1] == "Linux")
{
  restorePath = '/home/jeon/Dropbox/GitHub/RankConsistency'
} else {
  restorePath = 'C:/Users/Jeon/Dropbox/GitHub/RankConsistency'
}

save(file = paste0(restorePath,
                   '/result/sreal_traninig_', i_1,"_", i_2), 
     list = c('sr_est.list',   'bt_est.list', 
              'gbt_est.list', 'sr1_est.list',
              'gbt_est.list2',
              'sim.num'))

