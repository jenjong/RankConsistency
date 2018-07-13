# preprocessing function to result table to wide format



# preprocessing function to compute n_jk and w_jk
QmatFun <- function(race_mat, num_vec, p = 43, sel_idx = 1:43)
{
  
  n_mat <- matrix(0, p, p)  ## n_mat_jk 
  w_mat <- matrix(0, p, p)  ## w_mat_jk
  i=1
  for (i in 1:nrow(race_mat))  ## nrow 
  {
    n_v <- num_vec[i] ## n_v 
    p_vec<- race_mat[i,1:n_v] 
    for (j in 1:(n_v-1))
    {
      for (k in (j+1):n_v)
      {
        if (p_vec[j] == p_vec[k]) next  ## ties
        w_mat[ p_vec[j], p_vec[k] ] <- w_mat[ p_vec[j], p_vec[k] ] + 1
        n_mat[ p_vec[j], p_vec[k] ] <- n_mat[ p_vec[j], p_vec[k] ] + 1
        n_mat[ p_vec[k], p_vec[j] ] <- n_mat[ p_vec[j], p_vec[k] ]  
      }
    }
  }
  
  sc_alarm <- FALSE
  if ( min(rowSums(n_mat) ) == 0 )  
  {
    sc_alarm <- TRUE
  }
  Qmat <- n_mat[sel_idx, sel_idx]
  w_mat <- w_mat[sel_idx, sel_idx]
  
  
  n = sum(Qmat) 
  
  Gmat_hat <- w_mat/Qmat 
  Gmat_hat[!is.finite(Gmat_hat)] = 0  
  Qpmat = Qmat
  for (j in 1:nrow(Qmat)) Qpmat[j,] = Qmat[j,]/n*2  
  
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
  colnames(Gmat_hat) = colnames(Qpmat) = 
    colnames(Qmat) = colnames(w_mat) = names(sel_idx)
  
  return( list(x=x, y=y, Qpmat=Qpmat, Gmat_hat=Gmat_hat, 
               n = n, sc_alarm = sc_alarm, Qmat = Qmat,
               Wmat = w_mat) )
}

# Fit the Bradley model
btFun<- function(Qmat_fit)
{
  Qpmat = Qmat_fit$Qpmat  
  Gmat_hat = Qmat_fit$Gmat_hat
  x = Qmat_fit$x
  y = Qmat_fit$y
  p = ncol(Qpmat)  
  wmat = Qpmat*Gmat_hat
  wmat = t(wmat)
  wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  ## w_jj 제거 
  # fit glmnet
  fit <- glmnet(x, y, family = 'binomial',
                intercept = FALSE, weights = wvec, lambda = 0,
                standardize = F, 
                thresh = 1e-09)
  est = c(fit$beta[,1],0) ## lambda_43 추가 
  naive_est <- est
  names(naive_est) = colnames(Qpmat)
  return( naive_est )
}

# Fit the generalized Bradley model with cvec[k] (thresholding parameter)
gbtFun <-function(Qmat_fit, cut_v=0, ctype = 'boost')
{
  Qmat = Qmat_fit$Qmat
  Qpmat = Qmat_fit$Qpmat
  Gmat_hat = Qmat_fit$Gmat_hat

  sc_list = list()
  p = ncol(Qpmat)

  idx = 1
  result = matrix(0, p*(p-1)/2, 4)
  for ( i1 in 1:(p-1))
  {
    for (i2 in (i1+1):p) 
    {
      Qpmat.c1 = Qpmat
      idx1 <- ( Qpmat.c1[i1,] <= cut_v )
      idx2 <- ( Qpmat.c1[i2,] <= cut_v )
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
      Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cut_v)  
      nvec1 = Qpmat.c1[i1,]
      nvec2 = Qpmat.c1[i2,]
      idx1 = which(nvec1 == 0 | nvec2 == 0)
      idx2  =  setdiff( idx1, c(i1, i2))
      
      nvec3 = (nvec1[-idx1]+nvec2[-idx1])/2
      Qpmat.c2[i1,-idx1] = Qpmat.c2[i2,-idx1] = nvec3
      
      if (length(idx2)>0)  Qpmat.c2[i1,idx2] =  Qpmat.c2[i2,idx2] = 0
      
      Qpmat.c2[,i1] <- Qpmat.c2[i1,]
      Qpmat.c2[,i2] <- Qpmat.c2[i2,]  ## Qpmat.c2 : symm matrix
      # extend the graph
      {
      #idx3 <- sort( union(intersect( setdiff(1:p, idx1), setdiff(1:p, idx2) ),  c(i1, i2)) )
      ## find V_jk(maximum connected set)                
      i1i2_adj_matrix = matrix(as.integer(Qpmat.c2>0) , p , p)  ## adjacency matrix
      i1i2_graph = graph_from_adjacency_matrix(i1i2_adj_matrix , 
                                               mode="undirected" , weighted=NULL) 
      ## make a graph
      i1i2_clusters = clusters(i1i2_graph)$mem ## clustering using adj matrix
      if (i1i2_clusters[i1] != i1i2_clusters[i2])
      {  
        ## i1과 i2가 다른 connected 되지 않은 경우
        #  cat('   k:',k-1,', ',i1,'and',i2, 'is not connected!!\n')
        idx = idx + 1
        next  
      } 
      ## idx3 : edge index set of V_jk
      idx3 = sort(which(i1i2_clusters %in% i1i2_clusters[i1])) 
      }
      
      # correction
      
      if (ctype == 'balance')
      {
        j1_mat = Qpmat.c2[c(i1,i2),]
        j1_vec = j1_mat[j1_mat!=0]
        if (length(j1_vec)>0)
        {
          j1_var = mean(j1_vec)
          j1_mat[j1_mat!=0] = j1_var
          Qpmat.c2[c(i1,i2),] = j1_mat
        }
      }

      if (ctype == 'boost')
      {
        a = max(Qpmat.c2[c(i1,i2),])
        Qpmat.c2[c(i1,i2),][Qpmat.c2[c(i1,i2),]!=0] = a
      }
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
      
      pweight = sum(Qmat[idx3, idx3])/sum(Qmat)
      result[idx, 4] = pweight
      idx = idx + 1
    }
  }

  # recover gbt_est from sc_list
  sc_list <- result
  gbt_est = NULL
  
  tmp<-sc_list
  tmp <-tmp[tmp[,1]!=0, 1:4]
  p_set <-unique(c(tmp[,1:2]))
  if (length(p_set) != p) 
  {
    gbt_est = NULL
    cat('gbt_est is NULL!')
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
  names(gbt_est) = colnames(Qpmat)
  return( list(sc_list = sc_list, gbt_est = gbt_est) )
}

gbtFun_recov = function(result, Qmat_fit, method = 'binomial')
{
  sc_list <- result
  Qmat = Qmat_fit$Qmat
  Qpmat = Qmat_fit$Qpmat
  Gmat_hat = Qmat_fit$Gmat_hat
  p = ncol(Qpmat)
  
  if (method == "count")
  {
    mat = matrix(0,p,p);
    for (i in 1:nrow(result))
    {
      j = result[i,1];
      k = result[i,2];
      w_jk = result[i,3];
      if (j == 0) 
      {
        mat[j,k] = 0.5;
        mat[k,j] = 0.5;
        next;
      }
      if (w_jk==1) mat[j,k] = 1  else mat[k,j] = 1;
    }
    gbt_est = apply(mat,1,sum);
    return( gbt_est = gbt_est )
  }
    
  
  gbt_est = NULL
  tmp<-sc_list
  tmp <-tmp[tmp[,1]!=0, 1:4]
  p_set <-unique(c(tmp[,1:2]))
  if (length(p_set) != p) 
  {
    gbt_est = NULL
    cat('gbt_est is NULL!')
    return( gbt_est )
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
  fit<-glmnet(x,y, weights = w, family = method, lambda = 0)
  gbt_est <- c(fit$beta[,1],0)
  names(gbt_est) = colnames(Qpmat)
  return(gbt_est)
}

#rank_v<- c(1,2,4,3,6,5)
dcgFun <- function(rank_v)
{
  pp = length(rank_v)
  v = 0
  d_i = c(1,-diff(1/log(1:pp+1),lag=1))
  #d_i = rep(1,pp)
  #if (pp>5) d_i[(1:5)] = 0
  p_i = cumsum(d_i)
  #p_i = 1/log(2) + 1
  p_tmp = 1:pp - (1:pp)[rank_v]
  p_ib = abs(p_i - p_i[rank_v]) /abs(p_tmp)
  p_ib[is.nan(p_ib)] = 1
  for(i in 1:(pp-1))
  {
    for (j in (i+1):pp)
    {
      v <- v + p_ib[i]*p_ib[j]*as.integer(rank_v[i] > rank_v[j])
    }
  }
  return(v)
}


# measure the performance 
# evalFun_1: the average of kendall's tau from each game
evalFun_1 <- function(rdata, est, sel_idx)
{
  
  race_mat <- as.matrix(rdata[,18:33])
  num_vec<- rdata$V1
  perform_v <- rep(NA, length(num_vec))
  tmp = names(est)
  i = 1
  for (i in 1:length(num_vec))
  {
    obs_cars <- race_mat[i,][1:num_vec[i]]
    
    obs_idx = obs_cars %in% sel_idx
    if ( sum(obs_idx) < 2 ) next
    obs_cars <- obs_cars[obs_idx]
    if (length(unique(obs_cars)) <2) next
    
    rankest = length(est[match(obs_cars, sel_idx)]) - 
      rank( est[match(obs_cars, sel_idx)])  + 1
    if (length(unique(rankest)) == 1) next
    perform_v[i] = cor(1:length(obs_cars), rankest, method = 'kendall')
  }
  mean(perform_v, na.rm = T)  
}

# evalFun_2: the average of kendall's tau from total game
evalFun_2 <- function(rdata, est, sel_idx)
{
  
  race_mat <- as.matrix(rdata[,18:33])
  num_vec<- rdata$V1
  perform_v <- rep(NA, length(num_vec))
  tmp = names(est)
  i = 1
  n = 0
  for (i in 1:length(num_vec))
  {
    obs_cars <- race_mat[i,][1:num_vec[i]]
    
    obs_idx = obs_cars %in% sel_idx
    if ( sum(obs_idx) < 2 ) next
    obs_cars <- obs_cars[obs_idx]
    if (length(unique(obs_cars)) <2) next
    
    rankest = length(est[match(obs_cars, sel_idx)]) - 
      rank( est[match(obs_cars, sel_idx)])  + 1
    if (length(unique(rankest)) == 1) next
    perform_v[i] = cov(1:length(obs_cars), rankest, method = 'kendall')
    n = n + length(rankest)*(length(rankest)-1)
  }
  sum(perform_v, na.rm = T)/n 
}

# evalFun_3
evalFun_3 <- function(Qmat_fit, est)
{
  Gmat_hat = Qmat_fit$Gmat_hat
  Qmat = Qmat_fit$Qmat
  for ( i in 1:(length(est)-1))
  {
    for (j in (i+1):length(est))
    {
      if ( est[i]-est[j] < 0) Gmat_hat[i,j] = NA 
      if ( est[i]-est[j] > 0) Gmat_hat[j,i] = NA
      if ( est[i]== est[j])
      {
        Gmat_hat[j,i] =  Gmat_hat[i,j] = 0.5
      }
    }
  }
  Gmat_hat[Qmat==0] = NA
  
  sum(Gmat_hat, na.rm = T)/sum(!is.na(Gmat_hat))
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

# evalFun_4
evalFun_4 <- function(Qmat_fit, est)
{
  Gmat_hat = Qmat_fit$Gmat_hat
  Qmat = Qmat_fit$Qmat
  Gmat_hat[Qmat == 0] = NA
  GG = rep(0,length(est))
  
  for ( i in 1:(length(est)))
  {
    for (j in 1:length(est))
    {
      if (is.na(Gmat_hat[i,j])) next
      
      if ( est[i]- est[j] > 0) 
      {
        GG[i] = GG[i] + Gmat_hat[i,j]
      }
      
      if ( est[i]- est[j] < 0) 
      {
        GG[i] = GG[i] + Gmat_hat[j,i]
      }
      if ( est[i]== est[j])
      {
        GG[i] = GG[i] + 0.5
      }
    }
  }
  
  v1 = GG
  v2 = apply(!is.na(Gmat_hat), 1, sum, na.rm = T)
  names(v1) = names(v2)  = colnames(Gmat_hat)
  list(v1 = v1, v2 = v2)
}



evalFun_4_pair = function(result, Qmat_fit)
{
  tmp1 = result
  tmp2<-cbind(result[,2], result[,1], 1 - result[,3], result[,4])
  
  tmp = rbind(tmp1, tmp2)
  
  Gmat_hat = Qmat_fit$Gmat_hat
  Qmat = Qmat_fit$Qmat
  Gmat_hat[Qmat == 0] = NA
  GG = rep(0,ncol(Qmat))
  
  k = 1
  for ( k in 1:nrow(tmp))
  {
    i = tmp[k,1] 
    j = tmp[k,2]
    if (is.na(Gmat_hat[i,j])) next
    if (tmp[k,3] == 1  )   GG[i] = GG[i] + Gmat_hat[i,j]
    if (tmp[k,3] == 0  )   GG[i] = GG[i] + Gmat_hat[j,i]
  }
  
  v1 = GG
  v2 = apply(!is.na(Gmat_hat), 1, sum, na.rm = T)
  names(v1) = names(v2)  = colnames(Gmat_hat)
  list(v1 = v1, v2 = v2)
}

# sr
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
  result[,4] = 1
  result
}

# note that this sr uses only paired comparison results
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
  result[,4] = 1
  result
}

balFun = function(obs_cars, bt_est, Qpmat)
{
  pp = length(obs_cars)
  est_coef = bt_est[obs_cars]
  v = 0
  i.num = 0
  for(i in 1:(pp-1))
  {
    for (j in (i+1):pp)
    {
      if (est_coef[i] == est_coef[j]) next
      bal.var = 1/Qpmat[obs_cars[i], obs_cars[j]]
      v = v + bal.var*as.integer(est_coef[i] < est_coef[j])
      i.num = i.num + 1
    }
  }
  return (v)
}

kenFun = function(obs_cars, bt_est)
{
  pp = length(obs_cars)
  est_coef = bt_est[obs_cars]
  v = 0
  i.num = 0
  for(i in 1:(pp-1))
  {
    for (j in (i+1):pp)
    {
      if (est_coef[i] == est_coef[j]) next
      if (obs_cars[i] == 35 | obs_cars[j] == 35) next
      v = v + as.integer(est_coef[i] < est_coef[j])
      i.num = i.num + 1
    }
  }
  if (i.num > 0) return (v/i.num) else return(0)
}

cv_gbtFun <- function(rdata, cvec,  sample_idx, kfold = 5)
{
  cv_k = 1
  cv_err_DCG<- cv_err_kendall <- NULL
  sid <- sample(1:kfold, length(sample_idx), replace = TRUE)
  for (cv_k in 1:kfold)
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
    tmp = gbt_eval(sc_list, race_mat_test, num_vec_test, cvec, 
                   return_list = FALSE)$tau_result_vec
    
    
    cv_err_kendall <- rbind(cv_err_kendall, tmp[1,])
    cv_err_DCG <- rbind(cv_err_DCG, tmp[2,])
  }
  return(list(cv_err_DCG = cv_err_DCG, cv_err_kendall = cv_err_kendall))
}  

