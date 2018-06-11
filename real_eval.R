##### evaluation function 
rm(list = ls())
setwd("C:/Users/jeon/Documents/GitHub/RankConsistency")
load("real_0421-sc.Rdata")
load("real_0421.Rdata")
source('./lib/real_lib.R')

# result_matrix_DCG = matrix(0,500, length(cvec_r)+1)
# for (seed_v in 1:inner_iter)
# {
#   cat("iteration::", seed_v, '\n')
#   seed_v_i = (file_idx -1)*inner_iter + seed_v
#   set.seed(seed_v_i)
#   sc_list = vector(mode ='list', length = max_k)
#   ## 논문에 나온대로 7:3으로 뽑음. 
#   sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.8)))  
#   # cross validation : 여기서 sample 다시 생성해야 함!
#   race_mat <- as.matrix(rdata[sample_idx,18:33])   ## train set의 각 게임당 선택 차종 
#   num_vec<- rdata$V1[sample_idx]  ## 각 게임마다 참여한 유저 수 
#   naive_est <- result_list$naive[[seed_v]]
#   
#   race_mat_test<- as.matrix(rdata[-sample_idx,18:33])
#   num_vec_test <- rdata$V1[-sample_idx]
#   ######## evaluate performances of standard BT estimator ####
#   naive_fit <- naive_eval(race_mat_test,num_vec_test,
#                           naive_est)
#   result_matrix_DCG[seed_v_i, 1] <- naive_fit$tau_result_vec 
#   ######## evaluate performances of the two estimator ####    
#   
#   gbt_fit <- result_list$gbt[[seed_v]]
#   tmpvec = rep(NA, length(cvec_r))
#   i = 1
#   for (i in 1:length(cvec_r))
#   {
#     gbt_est <-gbt_fit$gbt_est_mat[i,]
#     if (any(is.na(gbt_est))) next
#     tmpvec[i] = naive_eval(race_mat_test,num_vec_test,
#                            gbt_est)$tau_result
#   }
#   result_matrix_DCG[seed_v_i, 2:(length(cvec)+1)]<- tmpvec
#   report_v <- colMeans(tau_result_matrix[1:seed_v,,drop = F], na.rm = T )
#   cat('now::::\n')
#   cat(round(report_v,5),'\n') 
# }


# select cross validation index
min_vec_DCG = c()
min_vec_kendall = c()
file_idx = 1
for (file_idx in 1:20)
{
  load( paste0('./result/real_data/real_cv_result_',file_idx,'.Rdata') )
  ii = 1
  for (ii in 1:inner_iter)
  {
    cv_fit <- cv_list[[ii]]
    cv_fit
    sel_i = (file_idx -1)*inner_iter + ii
    min_vec_DCG[sel_i] = which.min(colMeans(cv_fit$cv_err_DCG, na.rm = T))
    min_vec_kendall[sel_i] = 
      which.min(colMeans(cv_fit$cv_err_kendall, na.rm = T))
  }
}
#write.csv(min_vec_kendall,"min_vec_kendall.csv")
# DCG
aa = c()
for ( i in 1:nrow(result_matrix_DCG))
  aa[i] = result_matrix_DCG[i, min_vec_DCG[i]+1]
result_matrix_DCG


boxplot(result_matrix_DCG[,1], 
        result_matrix_DCG_sc[,3],
        aa, names = c("BT", "SC", "gBT"), 
        ylab = "generalized rank distance", col='lightblue')

t.test(result_matrix_DCG[,1]-aa)
aa = c()
for ( i in 1:nrow(result_matrix_kendall))
  aa[i] = result_matrix_kendall[i, min_vec_kendall[i]+1]

boxplot(result_matrix_kendall[,1], result_matrix_kendall_sc[,3], aa,
        names = c("BT", "SC", "gBT"), 
        ylab = "Kendall's distance", col='lightblue')




