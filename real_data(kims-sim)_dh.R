rm(list = ls())
gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
library(igraph)
library(MASS)
source('car_lib.r')
source('lib_rank.r')
require('glmnet')
source('sim.R')
rdata<-read.csv('racing_data.csv', header=F)

max.k = 3
file.idx = 1
inner.iter = 10
seed.value = 1

tau.result.matrix <- matrix(0, inner.iter, max.k + 2)

seed.value = 1
for ( seed.value in 1:inner.iter)
{
  seed.value.i = (file.idx -1)*inner.iter + seed.value
  set.seed(seed.value.i)
  sc.list = vector(mode ='list', length = max.k)

  sample.idx <- sort( sample(1:nrow(rdata), trunc(nrow(rdata)*0.7)))  ## 논문에 나온대로 7:3으로 뽑음. 
  
  # cross validation : 여기서 sample 다시 생성해야 함!
  r1.mat<- as.matrix(rdata[sample.idx,18:33])   ## train set의 각 게임당 선택 차종 
  num.vec<- rdata$V1[sample.idx]  ## 각 게임마다 참여한 유저 수 
  n.mat <- matrix(0, 43, 43)  ## n.mat_jk : j,k 차종의 비교 수(symm mat) 
  w.mat <- matrix(0, 43, 43)  ## w.mat_jk : j,k 차종의 승리 수 
  i=1
  for ( i in 1:nrow(r1.mat))  ## nrow : 게임 수 
  {
    n.v <- num.vec[i] ## n.v : i번째 게임의 참여 유저 수 
    p.vec<- r1.mat[i,1:n.v] ## i번째 게임의 차종(순위대로..) 
    for (j in 1:(n.v-1))
    {
      for (k in (j+1):n.v)
      {
        if (p.vec[j] == p.vec[k]) next  ## 같은 차종이 있을 경우는 제외. 
        w.mat[ p.vec[j], p.vec[k] ] <- w.mat[ p.vec[j], p.vec[k] ] + 1
        n.mat[ p.vec[j], p.vec[k] ] <- n.mat[ p.vec[j], p.vec[k] ] + 1
        n.mat[ p.vec[k], p.vec[j] ] <- n.mat[ p.vec[j], p.vec[k] ]  ## n.mat은 symm matrix 
      }
    }
  }
  ##########
  if ( min(rowSums(n.mat) ) == 0 )  ## 한 번도 경기를 한 적 없는 차종이 존재하는 경우.. 
  {
    sc.list = NULL
#    save.image(file = paste("C:/Users/uos_stat/Dropbox/A rank-consistency/prog/result/real-",
#                            seed.value,'.rdata', sep =''))  
    next
  }
  Qmat <- n.mat # n_{jk}
  n = sum(Qmat) ## total 비교 수의 2배 
  nj = colSums(Qmat)  ## 각 차종당 비교 수 
  #Qmat <- n.mat ## 위와 같은 명령이므로 제거 
  #Qmat.c <- Qmat/sum(Qmat)
  Qmat.c = Qmat/n ## Qmat.c_jk = n_jk / n
  Gmat.hat <- w.mat.r <- w.mat
  Gmat.hat <- Gmat.hat/Qmat ## Gmat.hat_jk = j가 k를 이긴 횟수 / n_jk 
  Gmat.hat[!is.finite(Gmat.hat)] = 0  ## j,k의 경기가 없을 경우 0으로.. 
  Qpmat = Qmat
  for (j in 1:nrow(Qmat)) Qpmat[j,] = Qmat[j,]/n*2  ## n은 total 비교 수의 2배이므로... 
  
  #############################################
  ## x matrix와 y vector 만들기 
  p = ncol(Qmat)
  x = matrix(0, p*(p-1), p)
  y = rep(0, p*(p-1) )
  ix = 1
  for (i in 1:p)
  {
    for (j in 1:p)
    {
      if (i == j) next
      jx1 = min(i,j)
      jx2 = max(i,j)
      x[ix,jx1] = 1; x[ix,jx2] = -1
      if (i<j) y[ix] = 1
      ix = ix + 1
    }
  }
  x = x[,-p]
#############################################

    Result = NULL
    Naive.list = list()
    Result.list = list()
    k = 0
    Qpmat.c1 = Qpmat  

    ##########################################################
    ######## naive BT fit
    wmat = Qpmat.c1*Gmat.hat
    wmat = t(wmat)
    wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  ## w_jj 제거 
    # fit glmnet
    fit <- glmnet(x, y, family = 'binomial',
                  intercept = FALSE, weights = wvec, lambda = 0, standardize = F, thresh = 1e-09)
    est = c(fit$beta[,1],0) ## lambda_43 추가 
    naive.est <- est

    ##########################################################
    ######## gBT fit
    
    ###############################
    # set weight-vector 
    #cvec <- unique(sort(Qmat.c))
    cvec <- (0:max.k)/n*2 ## cvec : threshold c vector
    
    for ( k in 0:max.k)
    {
      cat('thershold::', k, '\n')
      i1 = 1 ; i2 = 2
      idx = 1
      result = matrix(0, p*(p-1)/2, 4)
      for ( i1 in 1:(p-1))
      {
        for (i2 in (i1+1):p) 
        {
          Qpmat.c1 = Qpmat
          idx1 <- ( Qpmat.c1[i1,] <= cvec[k+1] )
          idx2 <- ( Qpmat.c1[i2,] <= cvec[k+1] )
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
          Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cvec[k+1])  

          nvec1 = Qpmat.c1[i1,]
          nvec2 = Qpmat.c1[i2,]
          idx1 = which(nvec1 == 0 | nvec2 == 0)
          idx2  =  setdiff( idx1, c(i1, i2))
          
          nvec3 = (nvec1[-idx1]+nvec2[-idx1])/2
          Qpmat.c2[i1,-idx1] = Qpmat.c2[i2,-idx1] = nvec3
          
          if (length(idx2)>0)  Qpmat.c2[i1,idx2] =  Qpmat.c2[i2,idx2] = 0
          
          Qpmat.c2[,i1] <- Qpmat.c2[i1,]
          Qpmat.c2[,i2] <- Qpmat.c2[i2,]  ## Qpmat.c2 : symm matrix
          
          #idx3 <- sort( union(intersect( setdiff(1:p, idx1), setdiff(1:p, idx2) ),  c(i1, i2)) )

          ## find V_jk(maximum connected set)                
          i1i2_adj_matrix = matrix(as.integer(Qpmat.c2>0) , p , p)  ## adjacency matrix
          i1i2_graph = graph_from_adjacency_matrix(i1i2_adj_matrix , mode="undirected" , weighted=NULL) ## make a graph
          i1i2_clusters = clusters(i1i2_graph)$mem ## clustering using adj matrix
          if (i1i2_clusters[i1] != i1i2_clusters[i2]){  ## i1과 i2가 다른 connected 되지 않은 경우
            cat('   k:',k,', ',i1,'and',i2, 'is not connected!!\n')
            idx = idx + 1
            next  
          } 
          ## idx3 : edge index set of V_jk
          idx3 = sort(which(i1i2_clusters %in% i1i2_clusters[i1])) 

          #########################################
          ## computing gBT estimator
          #########################################  
          wmat <- Qpmat.c2[idx3,idx3]*Gmat.hat[idx3, idx3]
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
                                       intercept = FALSE, weights = wvec, lambda = 0.0001, alpha = 0, standardize = F, thresh = 1e-09),
                         silent = T)
          if (class(try.fit)[1] == 'try-error')  
          {
            idx = idx + 1
            next
          }
                    
          est = c(fit$beta[,1],0)
          result[idx, 1:2] = c(i1, i2)
          if( est[which(idx3==i1)] > est[which(idx3==i2)]) result[idx, 3] = 1
                   
          idx = idx + 1
        }
      }
      sc.list[[k+1]] <- result
    }
    
    ####################################################################################
    ## 여기까지 읽음
    ## 5 and 40 is not connected!!
    ## 9 and 40 is not connected!!
    ## 14 and 40 is not connected!!
    ## 16 and 40 is not connected!!
    ## 21 and 40 is not connected!!
    ## 26 and 40 is not connected!!
    ## 27 and 40 is not connected!!
    ## 28 and 40 is not connected!!
    ## 30 and 40 is not connected!!
    ## 31 and 40 is not connected!!
    ## 36 and 40 is not connected!!
    ## 37 and 40 is not connected!!
    ## 40 and 42 is not connected!!
    ####################################################################################

    
    # save.image(file = paste("C:/Users/uos_stat/Dropbox/A rank-consistency/prog/result/real-",
    #                         seed.value,'.rdata', sep =''))
    # 
    
    ##### end of pairwise learning ######
    
    ### make the test set #####
    ## test set의 각 게임당 선택 차종 
    r2.mat<- as.matrix(rdata[-sample.idx,18:33])
    num.vec<- rdata$V1[-sample.idx]
    n.mat <- matrix(0, 43, 43)
    w.mat <- matrix(0, 43, 43)
    i=1
    for ( i in 1:nrow(r2.mat))
    {
      n.v <- num.vec[i]
      p.vec<- r2.mat[i,1:n.v]
      for (j in 1:(n.v-1))
      {
        for (k in (j+1):n.v)
        {
          if (p.vec[j] == p.vec[k]) next
          w.mat[ p.vec[j], p.vec[k] ] <- w.mat[ p.vec[j], p.vec[k] ] + 1
          n.mat[ p.vec[j], p.vec[k] ] <- n.mat[ p.vec[j], p.vec[k] ] + 1
          n.mat[ p.vec[k], p.vec[j] ] <- n.mat[ p.vec[j], p.vec[k] ]
        }
      }
    }
    
    ######## evaluate performances of standard BT estimator ####    
    
    naive.rankest <- 44 - rank(naive.est)
    perform1 <- rep(0, length(num.vec))
    for (i in 1:length(num.vec))
    {
      obs_cars <- r2.mat[i,][1:num.vec[i]]
      rank_true <- 1:length(obs_cars)
      rank_hat  <- order( naive.est[obs_cars], decreasing = T)
      perform1[i] <- cor(rank_true, rank_hat, method = "kendall")
    }
    tau.result.matrix[seed.value, 1] <- mean(perform1)
    
    # naive.tau <- 0
    # for ( i in 1:43)
    # {
    #   for (j in 1:43)
    #   {
    #     if ( n.mat[i,j] == 0 ) next ## 대결해본 적 없는 쌍은 제외 
    #     naive.tau <- naive.tau - sign(naive.rankest[i] - naive.rankest[j])*w.mat[i,j]
    #   }
    # }
    #tau.result.matrix[seed.value, 1] <- naive.tau/sum(n.mat)

    ######## evaluate performances of the two estimator ####    
    k = 0
    for (k in 0:max.k)
    {
      tmp<-sc.list[[k+1]]
      tmp <-tmp[tmp[,1]!=0, 1:3]
      p.set <-unique(c(tmp[,1:2]))
      if (length(p.set) != 43) 
      {
        tau.result.matrix[seed.value, k+1] <- NA
        next
      }
      
      x <- matrix(0, nrow(tmp)*2, 43)
      y <- rep(0, nrow(tmp)*2)
      i = 1
      for ( i in 1:nrow(tmp))
      {
        vec1<-tmp[i,1:2]; vec2<- tmp[i,3]
        x[2*(i-1)+1, vec1] <- c(1,-1) ; y[2*(i-1)+1] <- vec2
        x[2*i, vec1] <- c(-1,1) ; y[2*i] <- abs(vec2 - 1)
      }
      
      x<- x[,-43]
      
      fit<-glmnet(x,y, family = 'binomial', lambda = 0.000001)
      gbt.est <- c(fit$beta[,1],0)
      gbt.rankest <- 44 - rank(gbt.est)

      perform2 <- rep(0, length(num.vec))
      for (i in 1:length(num.vec))
      {
        obs_cars <- r2.mat[i,][1:num.vec[i]]
        rank_true <- 1:length(obs_cars)
        rank_hat  <- order( gbt.est[obs_cars], decreasing = T)
        perform2[i] <- cor(rank_true, rank_hat, method = "kendall")
      }
      tau.result.matrix[seed.value, k+2] <- mean(perform2)
      # gbt.tau <- 0 
      # #i = 24 ; j = 3
      # for ( i in 1:43)
      # {
      #   for (j in 1:43)
      #   {
      #     if ( n.mat[i,j] == 0 ) next
      #     gbt.tau <- gbt.tau - sign(gbt.rankest[i] - gbt.rankest[j])*w.mat[i,j]    
      #   }
      # }
      # tau.result.matrix[seed.value, k+1] <- gbt.tau/sum(n.mat)
    }
    report.v <- colMeans(tau.result.matrix[1:seed.value,,drop = F], na.rm = T )
    cat('now::::\n')
    cat(round(report.v,5),'\n')
}

#####################################################################################
a.mat <- tau.result.matrix[,1:30]
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

idx <-!is.na(tau.result.matrix[,20])

mean(a)
mean(a.mat[,1])

boxplot(tau.result.matrix[idx,1:20])
abline( h = median(tau.result.matrix[idx,1]))

#seed.value = 201
#load(file = paste("C:/Users/uos_stat/Dropbox/A rank-consistency/prog/result/real-",
#                  seed.value,'.rdata', sep =''))


boxplot(tau.result.matrix[,1:20])
