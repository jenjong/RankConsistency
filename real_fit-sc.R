rm(list = ls())
gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
#setwd("C:/Users/uos_stat/Documents/GitHub/RankConsistency")
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')

rdata<-read.csv('racing_data.csv', header=F)
max_k = 4
cvec_r <- seq(0, max_k, by = 2)
file_idx = 1
inner_iter = 500
seed_v = 1
result_matrix_kendall = matrix(0,inner_iter, length(cvec_r)+1)
result_matrix_DCG = matrix(0,inner_iter, length(cvec_r)+1)

result_list = list()
result_list$naive = vector(mode = 'list', length = inner_iter)
result_list$gbt = vector(mode = 'list', length = inner_iter)
result_list$sr = vector(mode = 'list', length = inner_iter)  

for ( seed_v in 1:inner_iter)
{
  cat("iteration::", seed_v, '\n')
  seed_v_i = (file_idx -1)*inner_iter + seed_v
  set.seed(seed_v_i)
  sc_list = vector(mode ='list', length = max_k)
  sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.8)))  
  # cross validation : 여기서 sample 다시 생성해야 함!
  race_mat <- as.matrix(rdata[sample_idx,18:33])   ## train set의 각 게임당 선택 차종 
  num_vec<- rdata$V1[sample_idx]  ## 각 게임마다 참여한 유저 수 
  Qmat_fit <-QmatFunc(race_mat, num_vec)  
  Qpmat = Qmat_fit$Qpmat  
  Gmat_hat = Qmat_fit$Gmat_hat
  Gmat_hat[Qpmat==0] = 0.5
  diag(Gmat_hat) = 0
  sc_est <- rowSums(Gmat_hat)
  result_list$sc_est[[seed_v]] <-  sc_est
  
  ### make the test set #####
  race_mat_test<- as.matrix(rdata[-sample_idx,18:33])
  num_vec_test <- rdata$V1[-sample_idx]
  ######## evaluate performances of standard BT estimator ####
  sc_fit <- naive_eval(race_mat_test,num_vec_test,
                          sc_est, return_list = FALSE)
  result_matrix_kendall[seed_v, 3] <- sc_fit$tau_result[1]
  result_matrix_DCG[seed_v, 3] <- sc_fit$tau_result[2]
}
#save.image("real_0421.rdata")
result_matrix_DCG_sc = result_matrix_DCG
result_matrix_kendall_sc = result_matrix_kendall
#save.image("real_0421-sc.rdata")

