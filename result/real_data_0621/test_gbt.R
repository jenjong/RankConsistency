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
max_k = 1
cvec_r <- seq(0, max_k, by = 2)
file_idx = 1
inner_iter = 1
seed_v = 1

#for ( seed_v in 1:inner_iter)
#{
  cat("iteration::", seed_v, '\n')
  seed_v_i = (file_idx -1)*inner_iter + seed_v
  set.seed(seed_v_i)
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
  cvec <- cvec_r/n*2 ## cvec : threshold c vector
  
  sc_listFun
  
  sc_list <- sc_listFun(cvec, Qpmat, Gmat_hat)
  
  
  
  sc_mat <-sc_list[[1]]
  
  n_mat = matrix(0,43,43)
  idx = 1
  for ( i in 1:42)
  {
    for (j in (i+1):43)
    {
      k1 = sc_mat[idx, 1]
      k2 = sc_mat[idx, 2]
      
      if (k1 == 0)
      {
        n_mat[i, j] = n_mat[j, i] = NA
        idx = idx + 1 
        next
      }
      
      if (sc_mat[idx, 3] == 1) 
      {
        n_mat[i, j] = 1
      } else {
        n_mat[j, i] = 1
      }
      
      idx = idx + 1    
    }
  
  }
  colSums(is.na(n_mat))
  44-rank(apply(n_mat,1,sum, na.rm = T))
  
  ### make the test set #####
  gbt_fit <- gbt_eval(sc_list, race_mat_test = NULL, num_vec_test = NULL, cvec, 
                      return_list = FALSE)
  
  
  plot(44-rank(gbt_fit$gbt_est_mat[1,]) , 44-rank(naive_est), pch = 19,
       col = heat.colors(43)[rank(apply(Qpmat,1,sum))]  )
  
  plot(44-rank(gbt_fit$gbt_est_mat[1,]),  44-rank(apply(n_mat,1,sum, na.rm = T)),pch = 19,
       col = heat.colors(43)[rank(apply(Qpmat,1,sum))]  )
  plot(44-rank(naive_est), 44-rank(apply(n_mat,1,sum, na.rm = T)))

    race_mat_test<- as.matrix(rdata[,18:33])
    num_vec_test <- rdata$V1
    Qmat_fit <-QmatFunc(race_mat_test, num_vec_test)  
    Qpmat = Qmat_fit$Qpmat  
    
    
    gbt_est  = gbt_fit$gbt_est_mat[1,]
    bt_est <- naive_est
    
    
    perform_kendall = matrix(0,length(num_vec_test),2)
    # evaluation
    i = 1
    for (i in 1:length(num_vec_test))
    {
      obs_cars <- race_mat_test[i,][1:num_vec_test[i]]
      rank_true <- 1:length(obs_cars)
      # BT
      rank_hat  <- order( bt_est[obs_cars], decreasing = T)
      #perform_kendall[i,1] <- (1-cor(rank_true, rank_hat, method = "kendall"))/2
      perform_kendall[i,1] <- kenFun(obs_cars, bt_est)
      # gBT
      rank_hat  <- order( gbt_est[obs_cars], decreasing = T)
      #perform_kendall[i,2] <- (1-cor(rank_true, rank_hat, method = "kendall"))/2
      perform_kendall[i,2] <- kenFun(obs_cars, gbt_est)
    }  
     apply(perform_kendall,2,mean, na.rm = TRUE)
     boxplot(perform_kendall)
    
  k2
  
  k1[,2]
  
  idx = k1[,1] > k1[,2]
  boxplot(k1[idx,c(2,1,3)], col = 'lightblue',
          names = c("gBT", "BT", "SC"), 
          ylab = "kendall's rank distance")
  
  boxplot(k3[,c(2,1,3)], col = 'lightblue',
          names = c("gBT", "BT", "SC"), 
          ylab = "generalized rank distance")
  
  
  seed_v = 1
  v1 = v2 = v3  = NULL
  for ( seed_v in 1:inner_iter)
  {
    gbt_est <-result_list$gbt[[seed_v]]$gbt_est_mat[min_vec[seed_v],]
    if ( any(is.na(gbt_est)) ) break
    bt_est <- result_list$naive[[seed_v]]
    sr_est <- result_list_sc$sc_est[[seed_v]]
    v1 = rbind(v1, 44-rank(bt_est) )
    v2 = rbind(v2, 44-rank(gbt_est) )
    v3 = rbind(v3, 44-rank(sr_est) )
  }
  
  idx = which(k1[,1] < k1[,2])
  
  
  boxplot(v1[idx,]-v2[idx,], col = topo.colors(100)[rank(colSums(Qpmat))],
          ylim = c(-12,9))
  abline(h = 0, col ='red')
  boxplot(v1[-idx,]-v2[-idx,], col = topo.colors(100)[rank(colSums(Qpmat))],
          ylim = c(-12,9))
  abline(h = 0, col ='red')
  abline(h = 0)
  which(k1[,1] < k1[,2])
  v1[3,] - v2[3,]
  v1_med = apply(v1,2, mean)
  v2_med = apply(v2,2, mean)
  plot(sort(v1_med))
  points(v2_med[order(v1_med)],col = "red")
  idx = order(v1_med)
  boxplot(v1[,idx], col  =  '#0000FF50')
  boxplot(v2[,idx], add = TRUE, col = '#FF000050', )
  abline(v = 31)
  abline(v = 4)
  idx[31]
  idx[1:10]
  sum(colSums(Qpmat)[idx[1:10]])/2
  sum(colSums(Qpmat))/2
  colSums(Qpmat)[35]/2
  which(idx==43)
  plot(colSums(Qpmat)/2)
  
  k1
  
  
  
#save.image("real_0614-1.rdata")




