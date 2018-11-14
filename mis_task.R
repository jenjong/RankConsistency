# mis task 1 (show the results of estimated ranks)
# my mistake
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
load('Real_BT_gBT2_cv5_all_data.rdata')
BT_est_rank
gBT2_est_rank
est_rank <- data.frame(BT_est_rank, gBT2_est_rank)
head(est_rank)
write.csv(est_rank, "est_rank.csv")


# mis task 2 (show the results of estimated ranks)
setwd("C:/Users/uos_stat/Documents/GitHub/RankConsistency")
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')
rdata<-read.csv('racing_data.csv', header=F)
race_mat <- as.matrix(rdata[,18:33])   ## train set의 각 게임당 선택 차종 
num_vec<- rdata$V1  ## 각 게임마다 참여한 유저 수 
Qmat_fit <-QmatFunc(race_mat, num_vec)  
Qmat = Qmat_fit$Qmat  
Wmat = Qmat_fit$Wmat
colnames(Qmat) = rownames(Qmat) = colnames(Wmat) = rownames(Wmat) =
  rownames(est_rank)
write.csv(Qmat, "Qmat.csv")
write.csv(Wmat, "Wmat.csv")

# mis task 3 (show the result based on cross matching)
rm(list = ls())
gc()
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
sim.num = 100

rdata<-read.csv('racing_data.csv', header=F)
n = nrow(rdata)

bt_est.list = gbt_est.list = sr1_est.list =
  sr_est.list = gbt_est.list2 =
  vector(mode='list', length = sim.num)

bt_result.list = gbt_result.list = sr1_result.list =
  sr_result.list = gbt_result.list2 =
  vector(mode='list', length = sim.num)
s_idx = 1:n
# training code
race_mat <- as.matrix(rdata[s_idx,18:33])
num_vec <- rdata$V1[s_idx]
Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                   p=43, sel_idx)
a = unname(BT_est_rank)
a = sort.int(a,index.return = T)$ix
Wmat = Qmat_fit$Wmat[a,a]
Qmat = Qmat_fit$Qmat[a,a]
r = Wmat/Qmat
r[is.nan(r)] = NA

r1 = r
tmp1 = Wmat
b1 = rowSums(Wmat)[1:13]
tmp1[1:13,1:13] = 0
b2 = rowSums(tmp1)[1:13]
b2/b1


# between
r1 = r
r1[1:13,1:13] = 0
l1 = apply(r1[,-(1:13)], 1, mean, na.rm=T)[1:13]
l1 = l1[!is.na(l1)]
round(l1,2)
cor(l1,12:1, method = 'kendall')

# write.csv(Qmat, "Qmat.csv")
# write.csv(Wmat, "Wmat.csv")
# write.csv(r, "prob.csv")


r2 = r
r2[1:13,-c(1:13)] = 0
l2 = apply(r2[,1:13], 1, mean, na.rm=T)[1:13]
l2 = l2[-9]

cor(l2, 12:1,method = 'kendall')
