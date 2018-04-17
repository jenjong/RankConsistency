rm(list = ls())
gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
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
inner_iter = 10
tau_result_matrix <- matrix(0, inner_iter, length(cvec_r)+1)

seed_v = 1

for ( seed_v in 1:inner_iter)
{
  seed_v_i = (file_idx -1)*inner_iter + seed_v
  set.seed(seed_v_i)
  sc_list = vector(mode ='list', length = max_k)
  sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.8)))  ## 논문에 나온대로 7:3으로 뽑음. 

    # cross validation : 여기서 sample 다시 생성해야 함!
  race_mat <- as.matrix(rdata[sample_idx,18:33])   ## train set의 각 게임당 선택 차종 
  num_vec<- rdata$V1[sample_idx]  ## 각 게임마다 참여한 유저 수 
  Qmat_fit <-QmatFunc(race_mat, num_vec)  
  Qpmat = Qmat_fit$Qpmat  
  Gmat_hat = Qmat_fit$Gmat_hat
  x = Qmat_fit$x
  y = Qmat_fit$y
  n = Qmat_fit$n
  ##########################################################
  ######## naive BT fit
  naive_est<- naive_btFunc(x,y, Qpmat, Gmat_hat)
  ##########################################################
  ######## gBT fit
  ###############################
  # set weight-vector 
  
  cvec <- cvec_r/n*2 ## cvec : threshold c vector
  sc_list <- sc_listFun(cvec, Qpmat, Gmat_hat)
  ##### end of pairwise learning ######
  ### make the test set #####
  ## test set의 각 게임당 선택 차종 
  race_mat_test<- as.matrix(rdata[-sample_idx,18:33])
  num_vec_test <- rdata$V1[-sample_idx]
  ######## evaluate performances of standard BT estimator ####    
  tau_result_matrix[seed_v, 1] <- naive_eval(race_mat_test,num_vec_test,
                                             naive_est)
######## evaluate performances of the two estimator ####    
  tau_result_matrix[seed_v, 2:(length(cvec)+1)]<- 
                gbt_eval(sc_list, race_mat_test, num_vec_test, cvec)
  report_v <- colMeans(tau_result_matrix[1:seed_v,,drop = F], na.rm = T )
  cat('now::::\n')
  cat(round(report_v,5),'\n')
}
    
### Cross validation

cvec_r <- seq(0, max_k, by = 5)
seed_v_i = (file_idx -1)*inner_iter + seed_v
set.seed(seed_v_i)
sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.7)))  ## 논문에 나온대로 7:3으로 뽑음. 
sid <- sample(1:5, length(sample_idx), replace = TRUE)
cv_k = 1
cv_err<- NULL
for (cv_k in 1:5)
{
  sample_idx_cvtr<- sample_idx[sid!=cv_k]
  race_mat <- as.matrix(rdata[sample_idx_cvtr,18:33])   ## train set의 각 게임당 선택 차종 
  num_vec<- rdata$V1[sample_idx_cvtr]  ## 각 게임마다 참여한 유저 수 
  Qmat_fit <-QmatFunc(race_mat, num_vec)  
  Qpmat = Qmat_fit$Qpmat  
  Gmat_hat = Qmat_fit$Gmat_hat
  x = Qmat_fit$x
  y = Qmat_fit$y
  n = Qmat_fit$n
  cvec <- cvec_r/n*2
  sc_list <- sc_listFun(cvec, Qpmat, Gmat_hat)
  ##### end of pairwise learning ######
  ### make the test set #####
  ## test set의 각 게임당 선택 차종 
  sample_idx_cvte<- sample_idx[sid==cv_k]
  race_mat_test<- as.matrix(rdata[sample_idx_cvte,18:33])
  num_vec_test <- rdata$V1[sample_idx_cvte]
  ######## evaluate performances of standard BT estimator ####    
  tmp = gbt_eval(sc_list, race_mat_test, num_vec_test, cvec)
  cv_err <- rbind(cv_err, tmp)
}

  
  


#####################################################################################
a.mat <- tau_result_matrix[,1:30]
a<-c()
a.idx.vec <- c()
i = 1
for( i in 1:200)
{
  if ( sum(is.na(a.mat[i,]))>0 ) a.idx <- min(which(is.na(a.mat[i,]))) - 1
  else a.idx <- ncol(a.mat)
  a.idx.vec[i] <- a.idx
  a[i] <- a.mat[i,a.idx]
}

boxplot(a.mat[,-1], col= 'lightblue', xlab = 'thresholding level', ylab = 'correlation')
abline( h= mean(a.mat[,1]), lty = 2)

plot(colMeans(!is.na(a.mat[,-1])), type = 'b', lty = 2, ylab = 'recovery probability',
     xlab = 'thresholding level')

boxplot(a.mat[,1], a, names = c('BT', 'gBT-BT'), ylab = 'correlation', col= 'lightblue')
barplot(table(a.idx.vec)/200, xlab = 'thresholding level', ylab = 'percentage')




t.test(a.mat[,1], a)

idx <-!is.na(tau_result_matrix[,20])

mean(a)
mean(a.mat[,1])

boxplot(tau_result_matrix[idx,1:20])
abline( h = median(tau_result_matrix[idx,1]))

#seed_v = 201
#load(file = paste("C:/Users/uos_stat/Dropbox/A rank-consistency/prog/result/real-",
#                  seed_v,'.rdata', sep =''))


boxplot(tau_result_matrix[,1:20])
