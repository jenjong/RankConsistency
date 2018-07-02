
rm(list = ls())
gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
#setwd("C:/Users/uos_stat/Documents/GitHub/RankConsistency")
load("Real_BT_gBT2_cv5_all_data.rdata")
which(BT_est_rank==9)
gBT2_est_rank[40]
which(BT_est_rank==4)
gBT2_est_rank[12]
sel_idx = which(BT_est_rank <=14)
if ( any(sel_idx==12 ) ) sel_idx = sel_idx[sel_idx!=12]
#sel_idx = sel_idx[-5]

library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')

rdata<-read.csv('racing_data.csv', header=F)
max_k = 0
cvec_r <- seq(0, max_k, by = 2)
file_idx = 1
inner_iter = 1
seed_v = 1
result_matrix_kendall = matrix(0,inner_iter, length(cvec_r)+1)
result_matrix_DCG = matrix(0,inner_iter, length(cvec_r)+1)

result_list = list()
result_list$naive = vector(mode = 'list', length = inner_iter)
result_list$gbt = vector(mode = 'list', length = inner_iter)
cat("iteration::", seed_v, '\n')
seed_v_i = (file_idx -1)*inner_iter + seed_v
set.seed(seed_v_i)
sc_list = vector(mode ='list', length = max_k)
## 논문에 나온대로 7:3으로 뽑음. 
sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*1)))  
# cross validation : 여기서 sample 다시 생성해야 함!
race_mat <- as.matrix(rdata[sample_idx,18:33])   ## train set의 각 게임당 선택 차종 
num_vec<- rdata$V1[sample_idx]  ## 각 게임마다 참여한 유저 수 
Qmat_fit <-QmatFunc(race_mat, num_vec, 43, sel_idx)  
Qpmat = Qmat_fit$Qpmat  
Gmat_hat = Qmat_fit$Gmat_hat

x = Qmat_fit$x
y = Qmat_fit$y
n = Qmat_fit$n
######## naive BT fit


naive_est <- naive_btFunc(x,y, Qpmat, Gmat_hat)
sc_list <- sc_listFun(0, Qpmat, Gmat_hat)
gbt_est <- gbt_recv(sc_list, p = length(sel_idx))
BT_rank_est = length(sel_idx) + 1 - rank(naive_est)
gbT2_rank_est = length(sel_idx) + 1 -rank(gbt_est)
names(gbT2_rank_est) <- names(sel_idx)
names(BT_rank_est) <- names(sel_idx)
plot(BT_rank_est, gbT2_rank_est)



Wmat = Qmat_fit$Wmat
Qmat = Qmat_fit$Qmat
Gmat = Qmat_fit$Gmat_hat
colnames(Wmat) =colnames(Qmat) = colnames(Gmat) = names(BT_rank_est)
cbind(BT_rank_est,gbT2_rank_est)

#ix = names(BT_rank_est)[order(gbT2_rank_est)]
ix = names(BT_rank_est)[order(BT_rank_est)]
i = 1
v = c()
ii= 1
for (i in 1:(length(ix)-1) )
{
  aa1 = which(names(BT_rank_est) == ix[i])
  Qmat[aa1,]
  Wmat[aa1,]
  Wmat[aa1,]/Qmat[aa1,]
  idx1 = Qmat[aa1,] != 0 
  
  aa2 = which(names(BT_rank_est) == ix[i+1])
  Qmat[aa2,]
  Wmat[aa2,]
  Wmat[aa2,]/Qmat[aa2,]
  idx2 = Qmat[aa2,] != 0 
  idx = idx1 & idx2

  v1 = sum(Wmat[aa1,idx]/Qmat[aa1,idx])
  v2 = sum(Wmat[aa2,idx]/Qmat[aa2,idx])
  v[i] = v1 - v2
}
sum(v>0)

## pairwise

aa1 = which(names(BT_rank_est) == 'car13')
Qmat[aa1,]
Wmat[aa1,]
Wmat[aa1,]/Qmat[aa1,]
idx1 = Qmat[aa1,] != 0 

aa2 = which(names(BT_rank_est) == 'car7')
Qmat[aa2,]
Wmat[aa2,]
Wmat[aa2,]/Qmat[aa2,]
idx2 = Qmat[aa2,] != 0 
idx = idx1 & idx2

sum(Wmat[aa1,idx]/Qmat[aa1,idx])
sum(Wmat[aa2,idx]/Qmat[aa2,idx])
cbind(BT_rank_est,gbT2_rank_est)

### sc_list



tmp<-sc_list[[1]]
tmp <-tmp[tmp[,1]!=0, 1:3]
p_set <-unique(c(tmp[,1:2]))
if (length(p_set) != length(sel_idx)) 
{
  tau_result[,k] <- NA
  next
}

p = length(sel_idx)
x <- matrix(0, nrow(tmp)*2, p)
y <- rep(0, nrow(tmp)*2)

for ( i in 1:nrow(tmp))
{
  vec1<-tmp[i,1:2]; vec2<- tmp[i,3]
  x[2*(i-1)+1, vec1] <- c(1,-1) ; y[2*(i-1)+1] <- vec2
  x[2*i, vec1] <- c(-1,1) ; y[2*i] <- abs(vec2 - 1)
}
x<- x[,-p]
fit<-glmnet(x, y, family = 'binomial', lambda = 0)
gbt_est <- c(fit$beta[,1],0)
gbt_est   
naive_est

# measure * 3 | 13 vs 43

