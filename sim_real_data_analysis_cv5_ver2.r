rm(list = ls())
gc()

library(igraph)
library(MASS)



################################################################################
################################################################################
## Real data simulation 2
################################################################################
################################################################################

setwd("E:\\rank consistency\\simulation\\")

source('.\\code\\car_lib.r')
source('.\\code\\lib_rank.r')
require('glmnet')
source('.\\code\\sim.R')
source('.\\code\\sim_functions.r')

## data 불러오기
rdata<-read.csv('.\\data\\racing_data.csv', header=F)

num_iter = 10

num_subjects=43
cvec = 7:11
n_folds = 5

BT_gBT2_matrix = matrix(0 , nrow=num_iter , ncol=2)

for (i in 1:num_iter){
  a = Sys.time()
  cat ('num of iteration :' , i , '\n')
  ## train, test sample 만들기
  set.seed(123+i+50)
  train_ind = sort(sample(1:nrow(rdata) , trunc(nrow(rdata)*0.7)))
  train_num_vec = rdata[train_ind , 1]
  train_data = as.matrix(rdata[train_ind , 18:33])
  test_num_vec = rdata[-train_ind , 1]
  test_data = as.matrix(rdata[-train_ind , 18:33])

  n = nrow(train_data)
  ## cv fold index set 만들기
  cv_index = get_fold_index(n , n_folds=n_folds , seed=(1234+i))
  
  cv_tau_matrix = matrix(0 , nrow=n_folds , ncol=(length(cvec)+1))
  for (j in 1:n_folds){
    cat ('    cv iteration :' , j , '\n')
    j_test_ind = cv_index[[j]]
    j_train_num_vec = train_num_vec[-j_test_ind]
    j_train_data = train_data[-j_test_ind,]
    j_test_num_vec = train_num_vec[j_test_ind]
    j_test_data = train_data[j_test_ind,]

    ## BT vs gBT
    j_tau_vec = get_BT_gBT2_tau(num_subjects , j_train_num_vec , j_train_data ,
                                j_test_num_vec , j_test_data , cvec)
                                
    ## 결과 저장
    cv_tau_matrix[j,] = j_tau_vec
  }
  
  ## optimal c 찾기
  #not_na_ind = rep(T , ncol(cv_tau_matrix))
  #for (k in 1:nrow(cv_tau_matrix)){
  #  k_not_na_ind = !is.na(cv_tau_matrix[k , ])
  #  not_na_ind = not_na_ind & k_not_na_ind
  #}
  #if (sum(not_na_ind)==0){
  #  cat ('      all rows have at least 1 NA\n')
  #  next
  #}
  #cv_mean_tau = colMeans(cv_tau_matrix[not_na_ind , ] , na.rm = T)
  cv_mean_tau = colMeans(cv_tau_matrix , na.rm = T)
  if (sum(is.nan(cv_mean_tau[-1]))==length(cvec)){
    cat ('      all rows have NA\n')
    next
  }
  optimal_c = cvec[which.max(cv_mean_tau[2:length(cv_mean_tau)])]
  cat ('      optimal threshold is' , optimal_c , '\n')
  
  ## optimal c를 이용하여 최종적으로 model을 fitting하고, 이를 이용하여 test set에서의 tau를 저장
  i_tau_vec = get_BT_gBT2_tau(num_subjects , train_num_vec , train_data ,
                              test_num_vec , test_data , optimal_c)
                              
  BT_gBT2_matrix[i,] = i_tau_vec
  b = Sys.time()
  
  cat('       Time required for iteration :' , b-a , '\n')
}

save.image('.\\Rdata\\Real_BT_gBT2_cv5_ver3_2-3.rdata')
write.table(BT_gBT2_matrix , file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_2-3.csv' , sep=',' ,
            row.names=F , col.names=F)


rm(list=ls())
gc()
################################################################################
## BT_gBT2_matrix들 불러오기
setwd("E:\\rank consistency\\simulation\\")

#BT_gBT2_matrix_1_1 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_1-1.csv' , sep=',' , header=F)
#BT_gBT2_matrix_1_2 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_1-2.csv' , sep=',' , header=F)
#BT_gBT2_matrix_2_1 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_2-1.csv' , sep=',' , header=F)
#BT_gBT2_matrix_2_2 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_2-2.csv' , sep=',' , header=F)
#BT_gBT2_matrix_3_1 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_3-1.csv' , sep=',' , header=F)
#BT_gBT2_matrix_3_2 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_3-2.csv' , sep=',' , header=F)


#BT_gBT2_matrix = rbind(BT_gBT2_matrix_1_1 , BT_gBT2_matrix_1_2 , BT_gBT2_matrix_2_1 ,
#                       BT_gBT2_matrix_2_2 , BT_gBT2_matrix_3_1 , BT_gBT2_matrix_3_2)

BT_gBT2_matrix_1_1 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_1-1.csv' , sep=',' , header=F)
BT_gBT2_matrix_1_2 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_1-2.csv' , sep=',' , header=F)
BT_gBT2_matrix_1_3 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_1-3.csv' , sep=',' , header=F)
BT_gBT2_matrix_2_1 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_2-1.csv' , sep=',' , header=F)
BT_gBT2_matrix_2_2 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_2-2.csv' , sep=',' , header=F)
BT_gBT2_matrix_2_3 = read.csv(file='.\\Resdata\\Real_BT_gBT2_cv5_ver3_2-3.csv' , sep=',' , header=F)


BT_gBT2_matrix = rbind(BT_gBT2_matrix_1_1 , BT_gBT2_matrix_1_2 , BT_gBT2_matrix_1_3 , 
                       BT_gBT2_matrix_2_1 , BT_gBT2_matrix_2_2 , BT_gBT2_matrix_2_3)

not_0_ind = which(BT_gBT2_matrix[,1]!=0)
length(not_0_ind)

colMeans(BT_gBT2_matrix[not_0_ind,])


## plot 그리기

tmp_matrix = as.data.frame(BT_gBT2_matrix[not_0_ind,])
names(tmp_matrix) = c('BT' , 'gBT')
setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim_real_data_2_ver3_dh.jpg')
boxplot(tmp_matrix, col= 'lightblue', ylab = 'rank-accuracy' , 
        cex.lab=2, cex.axis=2, cex.sub=2 , lwd=1.5)
dev.off()



################################################################################
################################################################################
## Real data simulation 3
################################################################################
################################################################################

rm(list = ls())
gc()

setwd("E:\\rank consistency\\simulation\\")

source('.\\code\\car_lib.r')
source('.\\code\\lib_rank.r')
require('glmnet')
source('.\\code\\sim.R')
source('.\\code\\sim_functions.r')

## data 불러오기
rdata<-read.csv('.\\data\\racing_data.csv', header=F)
num_subjects=43
cvec = 0:11
n_folds=5

train_num_vec = rdata[ , 1]
train_data = as.matrix(rdata[ , 18:33])


n = nrow(train_data)
## cv fold index set 만들기
cv_index = get_fold_index(n , n_folds=n_folds , seed=(1234))

cv_tau_matrix = matrix(0 , nrow=n_folds , ncol=(length(cvec)+1))

for (j in 1:n_folds){
  cat ('cv iteration :' , j , '\n')
  j_test_ind = cv_index[[j]]
  j_train_num_vec = train_num_vec[-j_test_ind]
  j_train_data = train_data[-j_test_ind,]
  j_test_num_vec = train_num_vec[j_test_ind]
  j_test_data = train_data[j_test_ind,]

  ## BT vs gBT
  j_results = get_BT_gBT2_results(num_subjects , j_train_num_vec , j_train_data ,
                              j_test_num_vec , j_test_data , cvec)
  j_tau_vec = j_results$tau_vec

  ## 결과 저장
  cv_tau_matrix[j,] = j_tau_vec
}

## optimal c 찾기
cv_mean_tau = colMeans(cv_tau_matrix , na.rm = T)
optimal_c = cvec[which.max(cv_mean_tau[2:length(cv_mean_tau)])]
cat ('  optimal threshold is' , optimal_c , '\n')

## optimal c를 이용하여 최종적으로 model을 fitting하고, 여기서 나오는 estimated rank를
## BT의 estimated rank와 비교
optimal_results = get_BT_gBT2_results(num_subjects , train_num_vec , train_data ,
                                      train_num_vec , train_data , optimal_c)

optimal_est = optimal_results$est_list

BT_est_rank = num_subjects+1-rank(optimal_est[[1]])
gBT2_est_rank = num_subjects+1-rank(optimal_est[[2]])

names(BT_est_rank) = names(gBT2_est_rank) = paste('car' , 1:43 , sep='')

BT_car_15 = names(sort(BT_est_rank)[1:15])
gBT2_car_15 = names(sort(gBT2_est_rank)[1:15])

match_ind = match(gBT2_car_15 , BT_car_15)

setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim_real_data_top15_ranks_dh.jpg')
plot(sort(gBT2_est_rank)[1:15] , match_ind , type='c' , lty=1 , xlab='gBT2 ranks' , 
     ylab='BT ranks' , xlim=c(0.6,15.5))
text(sort(gBT2_est_rank)[1:15] , match_ind , labels=gBT2_car_15 , col='red' , cex=1.3)
dev.off()

save.image('.\\Rdata\\Real_BT_gBT2_cv5_all_data.rdata')




################################################################################
## top 15 -> top 11
## 글씨 x , 점으로 표시, 선 x
## y=x 표시 
BT_car_11 = names(sort(BT_est_rank)[1:11])
gBT2_car_11 = names(sort(gBT2_est_rank)[1:11])

match_ind = match(gBT2_car_11 , BT_car_11)

setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim_real_data_top11_ranks_dh.jpg')
plot(sort(gBT2_est_rank)[1:11] , match_ind , type='p' , xlab='gBT ranks' , cex=2 , pch = 20 , 
     ylab='BT ranks' , xlim=c(0.6,11.5) , cex.lab=1.5, cex.axis=1.5, cex.sub=1.5)
abline(a=0,b=1 , cex.lab=2, cex.axis=2, cex.sub=2 , lwd=1.5)
dev.off()

