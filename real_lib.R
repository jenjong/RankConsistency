naive_btFunc<- function(x,y,Qpmat, Gmat_hat)
{
  p = ncol(Qpmat)  
  wmat = Qpmat*Gmat_hat
  wmat = t(wmat)
  wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  ## w_jj 제거 
  # fit glmnet
  fit <- glmnet(x, y, family = 'binomial',
                intercept = FALSE, weights = wvec, 
                lambda = 0, standardize = F, 
                thresh = 1e-09)
  est = c(fit$beta[,1],0) ## lambda_43 추가 
  naive_est <- est
  return( naive_est )
}



QmatFunc <- function(race_mat, num_vec)
{
  n_mat <- matrix(0, 43, 43)  ## n_mat_jk : j,k 차종의 비교 수(symm mat) 
  w_mat <- matrix(0, 43, 43)  ## w_mat_jk : j,k 차종의 승리 수 
  i=1
  for (i in 1:nrow(race_mat))  ## nrow : 게임 수 
  {
    n_v <- num_vec[i] ## n_v : i번째 게임의 참여 유저 수 
    p_vec<- race_mat[i,1:n_v] ## i번째 게임의 차종(순위대로..) 
    for (j in 1:(n_v-1))
    {
      for (k in (j+1):n_v)
      {
        if (p_vec[j] == p_vec[k]) next  ## 같은 차종이 있을 경우는 제외. 
        w_mat[ p_vec[j], p_vec[k] ] <- w_mat[ p_vec[j], p_vec[k] ] + 1
        n_mat[ p_vec[j], p_vec[k] ] <- n_mat[ p_vec[j], p_vec[k] ] + 1
        n_mat[ p_vec[k], p_vec[j] ] <- n_mat[ p_vec[j], p_vec[k] ]  ## n_mat은 symm matrix 
      }
    }
  }
  ##########
  sc_alarm <- FALSE
  if ( min(rowSums(n_mat) ) == 0 )  ## 한 번도 경기를 한 적 없는 차종이 존재하는 경우.. 
  {
    sc_alarm <- TRUE
  }
  Qmat <- n_mat # n_{jk}
  n = sum(Qmat) ## total 비교 수의 2배 
  
  Gmat_hat <- w_mat/Qmat ## Gmat_hat_jk = j가 k를 이긴 횟수 / n_jk 
  Gmat_hat[!is.finite(Gmat_hat)] = 0  ## j,k의 경기가 없을 경우 0으로.. 
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
  return( list(x=x, y=y, Qpmat=Qpmat, Gmat_hat=Gmat_hat, 
               n = n, sc_alarm = sc_alarm) )
}


sc_listFun<-function(cvec, Qpmat, Gmat_hat)
{
  sc_list = list()
  p = ncol(Qpmat)
  for ( k in 1:length(cvec))
  {
    #cat('thershold::', k-1, '\n')
    i1 = 1 ; i2 = 2
    idx = 1
    result = matrix(0, p*(p-1)/2, 4)
    for ( i1 in 1:(p-1))
    {
      for (i2 in (i1+1):p) 
      {
        Qpmat.c1 = Qpmat
        idx1 <- ( Qpmat.c1[i1,] <= cvec[k] )
        idx2 <- ( Qpmat.c1[i2,] <= cvec[k] )
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
        Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cvec[k])  
        
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
        i1i2_graph = graph_from_adjacency_matrix(i1i2_adj_matrix , 
                                                 mode="undirected" , weighted=NULL) ## make a graph
        i1i2_clusters = clusters(i1i2_graph)$mem ## clustering using adj matrix
        if (i1i2_clusters[i1] != i1i2_clusters[i2]){  ## i1과 i2가 다른 connected 되지 않은 경우
        #  cat('   k:',k-1,', ',i1,'and',i2, 'is not connected!!\n')
          idx = idx + 1
          next  
        } 
        ## idx3 : edge index set of V_jk
        idx3 = sort(which(i1i2_clusters %in% i1i2_clusters[i1])) 
        
        #########################################
        ## computing gBT estimator
        #########################################  
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
                                     lambda = 0.0001, alpha = 0, standardize = F,
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
        
        idx = idx + 1
      }
    }
    sc_list[[k]] <- result
  }
  return( sc_list )
}

naive_eval <- function(race_mat_test, num_vec_test, naive_est)
{
  naive_rankest <- 44 - rank(naive_est)
  perform_v <- rep(0, length(num_vec_test))
  for (i in 1:length(num_vec_test))
  {
    obs_cars <- race_mat_test[i,][1:num_vec_test[i]]
    rank_true <- 1:length(obs_cars)
    rank_hat  <- order( naive_est[obs_cars], decreasing = T)
    perform_v[i] <- cor(rank_true, rank_hat, method = "kendall")
  }
  return(mean(perform_v))
}

gbt_eval <- function(sc_list,race_mat_test, num_vec_test, cvec)
{
  tau_result_vec <- rep(0, length(cvec))
  gbt_est_mat <- matrix(NA, length(cvec), 43)
  for (k in 1:length(cvec))
  {
    tmp<-sc_list[[k]]
    tmp <-tmp[tmp[,1]!=0, 1:3]
    p_set <-unique(c(tmp[,1:2]))
    if (length(p_set) != 43) 
    {
      tau_result_vec[k] <- NA
      next
    }
    
    x <- matrix(0, nrow(tmp)*2, 43)
    y <- rep(0, nrow(tmp)*2)
    for ( i in 1:nrow(tmp))
    {
      vec1<-tmp[i,1:2]; vec2<- tmp[i,3]
      x[2*(i-1)+1, vec1] <- c(1,-1) ; y[2*(i-1)+1] <- vec2
      x[2*i, vec1] <- c(-1,1) ; y[2*i] <- abs(vec2 - 1)
    }
    x<- x[,-43]
    
    fit<-glmnet(x,y, family = 'binomial', lambda = 0.000001)
    gbt_est <- c(fit$beta[,1],0)
    gbt_est_mat[k,] <- gbt_est
    gbt_rankest <- 44 - rank(gbt_est)
    perform_v <- rep(0, length(num_vec_test))
    for (i in 1:length(num_vec_test))
    {
      obs_cars <- race_mat_test[i,][1:num_vec_test[i]]
      rank_true <- 1:length(obs_cars)
      rank_hat  <- order( gbt_est[obs_cars], decreasing = T)
      perform_v[i] <- cor(rank_true, rank_hat, method = "kendall")
    }
    tau_result_vec[k] <- mean(perform_v)
  }
  return(list(tau_result_vec = tau_result_vec, gbt_est_mat = gbt_est_mat))
}
