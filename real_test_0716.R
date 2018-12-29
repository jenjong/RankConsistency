rm(list = ls())
gc()
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/GitHub/RankConsistency")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
}

if (Sys.info()[1] == "Linux")
{
  restorePath = '/home/jeon/Dropbox/GitHub/RankConsistency'
} else {
  restorePath = 'C:/Users/Jeon/Dropbox/GitHub/RankConsistency'
}
# choose the training set
i_1 = 1
i_2 = 43
load(paste0(restorePath,
            '/result/real_traninig_', i_1,"_", i_2))
bt_est.list_1_43 = bt_est.list
gbt_est.list_1_43 = gbt_est.list
gbt_est.list2_1_43 = gbt_est.list2
sr1_est.list_1_43 = sr1_est.list
sr_est.list_1_43 = sr_est.list
# choose the top-tear 1-13
load("Real_BT_gBT2_cv5_all_data.rdata")
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
sel_idx.top = which(BT_est_rank >= 1 & BT_est_rank <= 13)
#plot(BT_est_rank,gBT2_est_rank, pch = 19, 
#     ylab = "ranks of gBT", xlab = "ranks of BT")
#ls()
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

if (Sys.info()[1] == "Linux")
{
  restorePath = '/home/jeon/Dropbox/GitHub/RankConsistency'
} else {
  restorePath = 'C:/Users/Jeon/Dropbox/GitHub/RankConsistency'
}
# overwritten
load(paste0(restorePath,
            '/result/sreal_traninig_', i_1,"_", i_2))
# test procedure
vmat1 = vmat2 = vmat3 = NULL
i = 1
#sel_idx.top = 1:43
sim.num = 50
for (i in 1:sim.num)
{
  cat(i,'\n')
  set.seed(i)
  s_idx = sample(1:n, trunc(n*0.7))
  
  race_mat <- as.matrix(rdata[-s_idx,18:33])
  num_vec <- rdata$V1[-s_idx]
  
  ### please note the use of sel_idx or sel_idx.top
  ### sel_idx.top is used when the full model is used.
  Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                     p=43, sel_idx.top)  
  #gbt_est = gbt_est.list[[i]]
  #if (is.null(gbt_est)) next
  gbt_est = gbt_est.list_1_43[[i]][sel_idx.top]
  #if (is.null(gbt_est)) next
  #bt_est = rank(bt_est.list[[i]])
  bt_est = bt_est.list_1_43[[i]][sel_idx.top]
  #sr_est = sr_est.list[[i]]
  sr_est = sr_est.list_1_43[[i]][sel_idx.top]
  #sr1_est = sr1_est.list[[i]]
  sr1_est = sr1_est.list_1_43[[i]][sel_idx.top]
  #gbt2_est = gbt_est.list[[i]]
  gbt2_est = gbt_est.list2_1_43[[i]][sel_idx.top]

  # v1 = evalFun_1(rdata[-s_idx,], bt_est, sel_idx)
  # v2 = evalFun_1(rdata[-s_idx,], gbt_est, sel_idx)
  # v5 = evalFun_1(rdata[-s_idx,], gbt2_est, sel_idx)
  # v3 = evalFun_1(rdata[-s_idx,], sr1_est, sel_idx)
  # v4 = evalFun_1(rdata[-s_idx,], sr_est, sel_idx)
  # vmat1 = rbind(vmat1, c(v1,v2,v5,v3,v4))
  # 
  # v1 = evalFun_2(rdata[-s_idx,], bt_est, sel_idx)
  # v2 = evalFun_2(rdata[-s_idx,], gbt_est, sel_idx)
  # v5 = evalFun_2(rdata[-s_idx,], gbt2_est, sel_idx)
  # v3 = evalFun_2(rdata[-s_idx,], sr1_est, sel_idx)
  # v4 = evalFun_2(rdata[-s_idx,], sr_est, sel_idx)
  # vmat2 = rbind(vmat2, c(v1,v2,v5,v3,v4))
  
  v1 = evalFun_3(Qmat_fit, bt_est)
  v2 = evalFun_3(Qmat_fit, gbt_est)
  v5 = evalFun_3(Qmat_fit, gbt2_est)
  v3 = evalFun_3(Qmat_fit, sr1_est)
  v4 = evalFun_3(Qmat_fit, sr_est)
  vmat3 = rbind(vmat3, c(v1,v2,v5,v3,v4))
}

# boxplot(vmat1[,-c(3)])
# boxplot(vmat2[,-c(3)])
boxplot(vmat3)
colMeans(vmat3)

if (Sys.info()[1] == "Linux")
{
  restorePath = '/home/jeon/Dropbox/GitHub/RankConsistency'
} else {
  restorePath = 'C:/Users/Jeon/Dropbox/GitHub/RankConsistency'
}

#save.image(file = paste0(restorePath,
#                         '/result/real_test_full'))


#save.image(file = paste0(restorePath,
#                   '/result/real_test_', i_1,"_", i_2))



rm(list = ls()) ; gc()

if (Sys.info()[1] == "Linux")
{
  restorePath = '/home/jeon/Dropbox/GitHub/RankConsistency'
} else {
  restorePath = 'C:/Users/Jeon/Dropbox/GitHub/RankConsistency'
}
# 7m 69
i_1 = 1; i_2 = 13
load(file = paste0(restorePath, '/result/real_test_', i_1,"_", i_2))
boxplot(vmat3[,c(5,1,4,2)], col='lightblue', 
        names = c('SC', 'BT', 'gSC', 'gBT'),
        ylab = 'accuracy')
colMeans(vmat3)[c(5,1,4,2)]





i_1 = 14; i_2 = 43
load(file = paste0(restorePath, '/result/real_test_', i_1,"_", i_2))
boxplot(vmat3[,c(5,1,4,2)], col='lightblue', 
        names = c('SC', 'BT', 'gSC', 'gBT'),
        ylab = 'accuracy')
colMeans(vmat3)

i_1 = 1; i_2 = 43
load(file = paste0(restorePath, '/result/real_test_', i_1,"_", i_2))
boxplot(vmat3[,c(4,1,3,2)], col='lightblue', 
        names = c('SC', 'BT', 'gSC', 'gBT'),
        ylab = 'accuracy')
colMeans(vmat3)
#
load(file = paste0(restorePath, '/result/real_test_full'))
boxplot(vmat3[,c(5,1,4,2)], col='lightblue', 
        names = c('SC', 'BT', 'gSC', 'gBT'),
        ylab = 'accuracy')
colMeans(vmat3)[c(5,1,4,2)]




