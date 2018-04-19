rm(list = ls())
gc()
#setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
setwd("C:/Users/uos_stat/Documents/GitHub/RankConsistency")
library(igraph)
library(MASS)
source('car_lib.R')
source('lib_rank.R')
source('sim.R')
source('real_lib.R')
require('glmnet')

rdata<-read.csv('racing_data.csv', header=F)
max_k = 15
cvec_r <- seq(0, max_k, by = 5)
file_idx = 1
inner_iter = 500
tau_result_matrix <- matrix(0, inner_iter, length(cvec_r)+1)

seed_v = 1
result_list = list()
result_list$naive = vector(mode = 'list', length = inner_iter)
result_list$gbt = vector(mode = 'list', length = inner_iter)
for ( seed_v in 1:inner_iter)
{
  cat("iteration::", seed_v, '\n')
  seed_v_i = (file_idx -1)*inner_iter + seed_v
  set.seed(seed_v_i)
  sc_list = vector(mode ='list', length = max_k)
  ## 논문에 나온대로 7:3으로 뽑음. 
  sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.8)))  
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
  ######## gBT fit
  
  cvec <- cvec_r/n*2 ## cvec : threshold c vector
  sc_list <- sc_listFun(cvec, Qpmat, Gmat_hat)
  ##### end of pairwise learning ######
  ### make the test set #####
  ## test set의 각 게임당 선택 차종 
  race_mat_test<- as.matrix(rdata[-sample_idx,18:33])
  num_vec_test <- rdata$V1[-sample_idx]
  ######## evaluate performances of standard BT estimator ####
  naive_fit <- naive_eval(race_mat_test,num_vec_test,
                             naive_est, return_list = TRUE)
  tau_result_matrix[seed_v, 1] <- naive_fit$tau_result_vec 
  
  ######## evaluate performances of the two estimator ####    
  gbt_fit <- gbt_eval(sc_list, race_mat_test, num_vec_test, cvec, 
                      return_list = TRUE)
  result_list$gbt[[seed_v]] <- gbt_fit
  tau_result_matrix[seed_v, 2:(length(cvec)+1)]<- gbt_fit$tau_result_vec
                
  report_v <- colMeans(tau_result_matrix[1:seed_v,,drop = F], na.rm = T )
  cat('now::::\n')
  cat(round(report_v,5),'\n')
}

# plot(44-rank(result_list$naive[[1]]), 44- rank(result_list$naive[[2]]))
# plot(44-rank(result_list$gbt[[1]]$gbt_est_mat[1,]), 
#              44- rank(result_list$naive[[2]]))

plot(naive_fit$perform_list + rnorm(1919,0,0.02), 
     gbt_fit$perform_list[[2]] + rnorm(1919,0,0.02),
     ylab = 'naive', xlab = 'gBT', col = "#00000010",
     pch = 20)

abline(a = 0 , b = 1)
str(gbt_fit$perform_list[[1]])
str(naive_fit$perform_list)
a = naive_fit$perform_list
b =  gbt_fit$perform_list[[1]]

race_mat_test[a==-1,]

race_mat_test[b==-1,]
