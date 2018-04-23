rm(list = ls())
gc()
#setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
setwd("~/task/RankConsistency")
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
file_idx = 3
inner_iter = 25
tau_result_matrix <- matrix(0, inner_iter, length(cvec_r)+1)

seed_v = 1
cv_list = vector(mode = 'list', length = inner_iter)
for ( seed_v in 1:inner_iter)
{
  cat("iteration::", seed_v, '\n')
  seed_v_i = (file_idx -1)*inner_iter + seed_v
  set.seed(seed_v_i)
  sc_list = vector(mode ='list', length = max_k)
  sample_idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.8)))  
  cv_fit <- cv_gbtFun(rdata, cvec,  sample_idx, kfold = 5)
  cv_list[[seed_v]] = cv_fit
}
save.image(paste0('real_cv_result_',file_idx,'.Rdata'))


