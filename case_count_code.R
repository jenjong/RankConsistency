############# Wd & library #############
rm(list = ls())
gc()

library(MASS)
library(igraph)
library(glmnet)
if (Sys.info()[1] == "Linux") {
  setwd("/home/jeon/Documents/GitHub/RankConsistency")
} else {
  setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
}
#setwd('C:\\Users\\moon\\Desktop\\racing_data')

# load car segmentation
load("Real_BT_gBT2_cv5_all_data.rdata")

source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')

sim.num = 100
rdata = read.csv('racing_data.csv', header=F)
n = nrow(rdata)

i_1 = 1; i_2 = 43;
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)

# training code
race_mat <- as.matrix(rdata[,18:33])
num_vec <- rdata$V1
Qmat_fit <-QmatFun(race_mat, num_vec, cut_var = 0,
                   p=43, sel_idx)  
Gmat_hat = Qmat_fit$Gmat_hat
Gmat_hat[Qmat_fit$Qmat==0] = NA
Gmat_hat[Gmat_hat==0.5] = NA


#### relation count
relation_count_fun = function(i, win_inx_mat, comb_mat)
{
  sub_win_inx_mat = win_inx_mat[comb_mat[i,], comb_mat[i,]]
  if(any(is.na(sub_win_inx_mat[row(sub_win_inx_mat) != col(sub_win_inx_mat)])))
  {
    return(NA) ## 관계 중 NA가 포함된 경우
  }else if(all(sort(rowSums(sub_win_inx_mat, na.rm = T)) == c(0, 1, 2))) ## 관계 성립 여부 확인
  {
    return(1) ## 관계 성립할 경우 1
  }else{
    return(0) ## 관계 성립하지 않을 경우 0
  }
}

relation_table_fun = function(Gmat_hat, inx_list)
{
  win_inx_mat = ifelse(Gmat_hat > 0.5, 1, 0)
  if(length(unique(lapply(inx_list, length))) == 1) ## case 1, 2, 3
  {
    comb_mat = t(combn(inx_list$i_inx, 3))
  }else{ ## case 4 - 5
    comb_expand_mat = expand.grid(inx_list$i_inx, inx_list$j_inx, inx_list$k_inx)
    comb_expand_list = apply(comb_expand_mat, 1, function(x) unique(sort(x)))
    if(class(comb_expand_list) == 'matrix')
    {
      comb_mat = t(comb_expand_list)
      comb_mat = comb_mat[!duplicated(comb_mat),]
    }else{
      comb_mat = comb_expand_list[unlist(lapply(comb_expand_list, length)) == 3] %>%
        do.call('rbind', .)
      comb_mat = comb_mat[!duplicated(comb_mat), ]
    }
  }
  relation_count_list = lapply(1:nrow(comb_mat), relation_count_fun, 
                               win_inx_mat = win_inx_mat, comb_mat = comb_mat)
  relation_count = unlist(relation_count_list)
  return(table(factor(relation_count, exclude = "")))
}

# Case 1
# 전체는 43등
i_1 = 1; i_2 = 43;
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
inx_list = list(
  i_inx = sel_idx,
  j_inx = sel_idx,
  k_inx = sel_idx
)
case1_result = relation_table_fun(Gmat_hat = Gmat_hat, inx_list = inx_list)
case1_result

# Case 2
# 1등에서 13등까지
i_1 = 1; i_2 = 13;
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
inx_list = list(
  i_inx = sel_idx,
  j_inx = sel_idx,
  k_inx = sel_idx
)
case2_result = relation_table_fun(Gmat_hat = Gmat_hat, inx_list = inx_list)
case2_result

# Case 3
# 14등에서 43등까지
i_1 = 1; i_2 = 13;
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
sel_idx_c = c(1:43)[-sel_idx]
inx_list = list(
  i_inx = sel_idx_c,
  j_inx = sel_idx_c,
  k_inx = sel_idx_c
)
case3_result = relation_table_fun(Gmat_hat = Gmat_hat, inx_list = inx_list)
case3_result

# Case 4
# i, j는 1등부터 13등까지 k는 14등부터 43등까지
i_1 = 1; i_2 = 13;
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
sel_idx_c = c(1:43)[-sel_idx]
inx_list = list(
  i_inx = sel_idx,
  j_inx = sel_idx,
  k_inx = sel_idx_c
)
case4_result = relation_table_fun(Gmat_hat = Gmat_hat, inx_list = inx_list)
case4_result

# Case 5
# i, j는 1등부터 13등까지 k는 14등부터 43등까지
i_1 = 1; i_2 = 13;
sel_idx = which(BT_est_rank >= i_1 & BT_est_rank <= i_2)
sel_idx_c = c(1:43)[-sel_idx]
inx_list = list(
  i_inx = sel_idx,
  j_inx = sel_idx_c,
  k_inx = sel_idx_c
)
case5_result = relation_table_fun(Gmat_hat = Gmat_hat, inx_list = inx_list)
case5_result

result_mat = rbind(case1_result,
                   case2_result,
                   case3_result,
                   case4_result,
                   case5_result)

write.csv(addmargins(result_mat, margin = 2), 'case_result_mat.csv')
