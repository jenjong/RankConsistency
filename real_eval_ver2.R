# 
# 1. load the list of coefficients from real_0421.rdata
# 2. there are bt and gbt estimators in result_list$naive and
#     result_list$gbt
# 3. 

# function: dcgFun()


rm(list = ls())
gc()
setwd("C:/Users/jeon/Documents/GitHub/RankConsistency")
load("real_0421-sc.Rdata")
result_list_sc = result_list
load("real_0421.rdata")

source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')

rdata<-read.csv('racing_data.csv', header=F)
min_vec_kendall<- read.csv("min_vec_kendall.csv")
min_vec<- min_vec_kendall$x
seed_v = 1
inner_iter = 100
file_idx = 1

k1 = matrix(0,inner_iter,3)
k2 = matrix(0,inner_iter,3)
k3 = matrix(0,inner_iter,3)

for ( seed_v in 1:inner_iter)
{
  cat("iteration::", seed_v, '\n')
  seed_v_i = (file_idx - 1)*inner_iter + seed_v
  set.seed(seed_v_i)
  sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.8)))  
  race_mat_test<- as.matrix(rdata[-sample_idx,18:33])
  num_vec_test <- rdata$V1[-sample_idx]
  Qmat_fit <-QmatFunc(race_mat_test, num_vec_test)  
  Qpmat = Qmat_fit$Qpmat  
  
  gbt_est <-result_list$gbt[[seed_v]]$gbt_est_mat[min_vec[seed_v],]
  bt_est <- result_list$naive[[seed_v]]
  sr_est <- result_list_sc$sc_est[[seed_v]]
  
  perform_kendall = matrix(0,length(num_vec_test),3)
  perform_gkendall = matrix(0,length(num_vec_test),3)
  perform_bkendall = matrix(0,length(num_vec_test),3)
  
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
    perform_gkendall[i,1] <-  dcgFun(rank_hat)
    perform_bkendall[i,1] <-  balFun(obs_cars, bt_est, Qpmat)
    
    # gBT
    rank_hat  <- order( gbt_est[obs_cars], decreasing = T)
    #perform_kendall[i,2] <- (1-cor(rank_true, rank_hat, method = "kendall"))/2
    perform_kendall[i,2] <- kenFun(obs_cars, gbt_est)
    perform_gkendall[i,2] <-  dcgFun(rank_hat)
    perform_bkendall[i,2] <-  balFun(obs_cars, gbt_est, Qpmat)
    # SR
    rank_hat  <- order( sr_est[obs_cars], decreasing = T)
    #perform_kendall[i,3] <- (1-cor(rank_true, rank_hat, method = "kendall"))/2
    perform_kendall[i,3] <- kenFun(obs_cars, sr_est)    
    perform_gkendall[i,3] <-  dcgFun(rank_hat)
    perform_bkendall[i,3] <-  balFun(obs_cars, sr_est, Qpmat)
  }  
  k1[seed_v_i,] = apply(perform_kendall,2,mean, na.rm = TRUE)
  k2[seed_v_i,] = apply(perform_gkendall,2, mean, na.rm = TRUE)
  k3[seed_v_i,] = apply(perform_bkendall,2, mean, na.rm = TRUE)
}  

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
