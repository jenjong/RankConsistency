


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
