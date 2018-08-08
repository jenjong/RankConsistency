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
