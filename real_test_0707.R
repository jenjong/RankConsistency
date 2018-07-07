rm(list = ls())
setwd('C:/Users/Jeon/Documents/GitHub/RankConsistency')
load("Real_BT_gBT2_cv5_all_data.rdata")
sel_idx = which(BT_est_rank <=13)
library(igraph)
library(MASS)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
require('glmnet')
rdata<-read.csv('racing_data.csv', header=F)
# data preprocessing
race_mat <- as.matrix(rdata[,18:33])
num_vec <- rdata$V1
Qmat_fit <-QmatFun(race_mat, num_vec, p=43, sel_idx) 
# pair 별로 추정이 잘 되는가?
sr1_fun = function(Qmat_fit)
{
  Gmat = Qmat_fit$Gmat_hat
  Qmat = Qmat_fit$Qmat
  Wmat = Qmat_fit$Wmat
  p = nrow(Gmat)
  result = matrix(0, p*(p-1)/2, 4)
  i = 1
  j = 8
  idx = 1
  for ( i in 1:(p-1))
  {
    for (j in (i+1):p)
    {
      n1 = Qmat[i,] ; n2 = Qmat[j,] ; v1 = Gmat[i,] ; v2 = Gmat[j,]
      withvec = n1*n2 != 0
      
      if ((Qmat[i,j] == 0) & (sum(withvec)==0)) 
      {
        idx = idx + 1
        next
      }
      
      if (Qmat[i,j] == 0) 
      {
        wr1 = mean(v1[withvec])
        wr2 = mean(v2[withvec])
      }
      
      if (Qmat[i,j] != 0) 
      {
        withvec[c(i,j)] = TRUE
        wr1 = sum(v1[withvec])/(length(v1[withvec])-1)
        wr2 = sum(v2[withvec])/(length(v2[withvec])-1)
      }
      
      result[idx, 1:3] = c(i,j, as.integer(wr1 > wr2) )
      idx = idx + 1
    }
  }
  result
}

# pair 
sr2_fun = function(Qmat_fit)
{
  Gmat = Qmat_fit$Gmat_hat
  Qmat = Qmat_fit$Qmat
  Wmat = Qmat_fit$Wmat
  p = nrow(Gmat) 
  result = matrix(0, p*(p-1)/2, 4)
  i = 1
  j = 2
  idx = 1
  for ( i in 1:(p-1))
  {
    for (j in (i+1):p)
    {
      n1 = Qmat[i,] ; n2 = Qmat[j,] ; v1 = Wmat[i,] ; v2 = Wmat[j,]
      wr1 = sum(v1)/sum(n1)
      wr2 = sum(v2)/sum(n2)
      result[idx, 1:3] = c(i,j, as.integer(wr1 > wr2) )
      idx = idx + 1
    }
  }
  result
}

evalFun_3_pair = function(result, Qmat_fit)
{
  Gmat_hat = Qmat_fit$Gmat_hat
  Qmat = Qmat_fit$Qmat
  p = ncol(Qmat)
  idx = 1
  for (i in 1:(p-1))
  {
    for ( j in (i+1):p)
    {
      if (result[idx,1] == 0) 
      {
        Gmat_hat[i,j] = Gmat_hat[j,i] = NA
        idx = idx + 1
        next
      }
      
      if ( result[idx,3] == 0) Gmat_hat[i,j] = NA else Gmat_hat[j,i] = NA
      idx = idx + 1
    }
  }
  Gmat_hat[Qmat==0] = NA
  sum(Gmat_hat, na.rm = T)/sum(!is.na(Gmat_hat))
}

###########################################################
result1 = sr1_fun(Qmat_fit)
result2 = sr2_fun(Qmat_fit)
gbt_fit <- gbtFun(Qmat_fit, cvec = 0)
result4 = gbt_fit$sc_list
evalFun_3_pair(result4, Qmat_fit)
  

# check recover
result1[,4] = 1
sc_list <- result1
gbt_est = NULL
tmp<-sc_list
tmp <-tmp[tmp[,1]!=0, 1:4]
p_set <-unique(c(tmp[,1:2]))
if (length(p_set) != p) 
{
  gbt_est = NULL
  cat('gbt_est is NULL!\n')
  return( list(sc_list = sc_list, gbt_est = gbt_est) )
}

x <- matrix(0, nrow(tmp)*2, p)
y <- rep(0, nrow(tmp)*2)
for ( i in 1:nrow(tmp))
{
  vec1<-tmp[i,1:2]; vec2<- tmp[i,3]
  x[2*(i-1)+1, vec1] <- c(1,-1) ; y[2*(i-1)+1] <- vec2
  x[2*i, vec1] <- c(-1,1) ; y[2*i] <- abs(vec2 - 1)
}
x<- x[,-p]
w = rep(tmp[,4], each = 2)
fit<-glmnet(x,y, weights = w, family = 'binomial', lambda = 0.000001)
gbt_est <- c(fit$beta[,1],0)
names(gbt_est) = colnames(Qmat)
gbt_est

evalFun_3(Qmat_fit, gbt_est)

## check the naive 
# compare the result of SR1 and gBT-pre
tmp = cbind(result1[,1:2], result1[,3], result4[,3])
idx = (result1[,3]!=result4[,3])
tmp[idx,]





# gbt code
cvec=0
Qmat = Qmat_fit$Qmat
Qpmat = Qmat_fit$Qpmat
Gmat_hat = Qmat_fit$Gmat_hat
sc_list = list()
p = ncol(Qpmat)
idx = 1
result = matrix(0, p*(p-1)/2, 4)
i1 = 2
i2 = 11

Qpmat.c1 = Qpmat
idx1 <- ( Qpmat.c1[i1,] <= cvec )
idx2 <- ( Qpmat.c1[i2,] <= cvec )
if (sum(idx1)>0 )
{
  Qpmat.c1[i1,idx1] <- 0 ;  Qpmat.c1[idx1,i1] <- 0
}
    
if (sum(idx2)>0 ) 
{
  Qpmat.c1[i2,idx2] <- 0 ;  Qpmat.c1[idx2,i2] <- 0  
}
    
Qpmat.c2 = Qpmat.c1
## thresholding procedure
Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cvec)  
nvec1 = Qpmat.c1[i1,]
nvec2 = Qpmat.c1[i2,]
idx1 = which(nvec1 == 0 | nvec2 == 0)
idx2  =  setdiff( idx1, c(i1, i2))
nvec3 = (nvec1[-idx1]+nvec2[-idx1])/2
Qpmat.c2[i1,-idx1] = Qpmat.c2[i2,-idx1] = nvec3
if (length(idx2)>0)  Qpmat.c2[i1,idx2] =  Qpmat.c2[i2,idx2] = 0
Qpmat.c2[,i1] <- Qpmat.c2[i1,]
Qpmat.c2[,i2] <- Qpmat.c2[i2,]  ## Qpmat.c2 : symm matrix
# extend graph

    #idx3 <- sort( union(intersect( setdiff(1:p, idx1), setdiff(1:p, idx2) ),  c(i1, i2)) )
    ## find V_jk(maximum connected set)                
    i1i2_adj_matrix = matrix(as.integer(Qpmat.c2>0) , p , p)  ## adjacency matrix
    i1i2_graph = graph_from_adjacency_matrix(i1i2_adj_matrix , 
                                             mode="undirected" , weighted=NULL) 
    ## make a graph
    i1i2_clusters = clusters(i1i2_graph)$mem ## clustering using adj matrix
    if (i1i2_clusters[i1] != i1i2_clusters[i2]){  
      ## i1과 i2가 다른 connected 되지 않은 경우
      #  cat('   k:',k-1,', ',i1,'and',i2, 'is not connected!!\n')
      idx = idx + 1
      next  
    } 
    ## idx3 : edge index set of V_jk
    idx3 = sort(which(i1i2_clusters %in% i1i2_clusters[i1])) 
    
    #########################################
    ## computing gBT estimator
    #########################################  
    j1_mat = Qpmat.c2[c(i1,i2),]
    if (j1_mat[1,i2]!=0)
    {
      a1 = j1_mat[1,i2]
      a2 = max(j1_mat[1,][j1_mat[1,]!=0])
      if (a1<a2) Qpmat.c2[i1,i2] = Qpmat.c2[i2,i1] = a2
    }
    wmat <- Qpmat.c2[idx3,idx3]*Gmat_hat[idx3, idx3]
    wmat = t(wmat)
    pp <- length(idx3)
    wvec = wmat[ - (1 + ( 0:(pp-1) ) *(pp+1))] 
    xx = matrix(0, pp*(pp-1), pp)
    yy = rep(0, pp*(pp-1) )
    ix = 1
    for (i in 1:pp)
    {
      for (j in 1:pp)
      {
        if (i == j) next
        jx1 = min(i,j)
        jx2 = max(i,j)
        xx[ix,jx1] = 1; xx[ix,jx2] = -1
        if (i<j) yy[ix] = 1
        ix = ix + 1
      }
    }
    xx = xx[,-pp]
    
    try.fit <- try(fit <- glmnet(xx, yy, family = 'binomial',
                                 intercept = FALSE, weights = wvec, 
                                 lambda = 0.0000, alpha = 0, standardize = F,
                                 thresh = 1e-09),
                   silent = T)
    if (class(try.fit)[1] == 'try-error')  
    {
      idx = idx + 1
      next
    }
    
    est = c(fit$beta[,1],0)
    result[idx, 1:2] = c(i1, i2)
    if( est[which(idx3==i1)] > est[which(idx3==i2)]) result[idx, 3] = 1
    
    pweight = sum(Qmat[idx3, idx3])/sum(Qmat)
    result[idx, 4] = pweight
    idx = idx + 1

