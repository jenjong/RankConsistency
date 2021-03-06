rm(list = ls())
gc()
if (Sys.info()[1] == "Linux" ) setwd("/home/jeon/Documents/Github/RankConsistency/result/real_data_0621")
if (Sys.info()[1] == "Windows" ) setwd("C:/Users/jeon/Documents/GitHub/RankConsistency/result/real_data_0621")
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')

rdata<-read.csv('racing_data.csv', header=F)
# argument data by 2 times
max_k = 10
cvec_r <- seq(0, max_k, by = 2)
file_idx = 1
inner_iter = 1
seed_v = 1

for ( seed_v in 1:inner_iter)
{
  cat("iteration::", seed_v, '\n')
  seed_v_i = (file_idx -1)*inner_iter + seed_v
  set.seed(seed_v_i)
  sc_list = vector(mode ='list', length = max_k)
  sample_idx <- 1:nrow(rdata)
  cv_fit <- cv_gbtFun(rdata, cvec_r,  sample_idx, kfold = 5)
  
  
  sc_list = vector(mode ='list', length = max_k)
  sample_idx <- 1:nrow(rdata)
  # cross validation : 여기서 sample 다시 생성해야 함!
  race_mat <- as.matrix(rdata[sample_idx,18:33])   ## train set의 각 게임당 선택 차종 
  num_vec<- rdata$V1[sample_idx]  ## 각 게임마다 참여한 유저 수 
  Qmat_fit <-QmatFunc(race_mat, num_vec)  
  Qpmat = Qmat_fit$Qpmat  
  Gmat_hat = Qmat_fit$Gmat_hat
  x = Qmat_fit$x
  y = Qmat_fit$y
  n = Qmat_fit$n
  ######## naive BT fit
  naive_est <- naive_btFunc(x,y, Qpmat, Gmat_hat)
  result_list$naive[[seed_v]] <-  naive_est
  cvec <- cvec_r/n*2 ## cvec : threshold c vector
  
  sc_list <- sc_listFun(cvec, Qpmat, Gmat_hat)

  ### make the test set #####
  race_mat_test<- as.matrix(rdata[-sample_idx,18:33])
  num_vec_test <- rdata$V1[-sample_idx]
  
  
  ######## evaluate performances of standard BT estimator ####
  naive_fit <- naive_eval(race_mat_test,num_vec_test,
                             naive_est, return_list = FALSE)
  result_matrix_kendall[seed_v, 1] <- naive_fit$tau_result[1]
  result_matrix_DCG[seed_v, 1] <- naive_fit$tau_result[2]
  ######## evaluate performances of the two estimator ####    
  gbt_fit <- gbt_eval(sc_list, race_mat_test = NULL, num_vec_test = NULL, cvec, 
                      return_list = FALSE)
  result_list$gbt[[seed_v]] <- gbt_fit
  result_matrix_kendall[seed_v, 2:(length(cvec)+1)]<- gbt_fit$tau_result_vec[1,]
  result_matrix_DCG[seed_v, 2:(length(cvec)+1)]<- gbt_fit$tau_result_vec[2,]
  report_v <- colMeans(result_matrix_DCG[1:seed_v,,drop = F], na.rm = T )
  cat('now::::\n')
  cat(round(report_v,5),'\n')
}

#save.image("real_0614-1.rdata")
list.files()

