#1  Dirichlet density function 
ddirichlet<- function(x,alpha.vec)
{
  xd = 1-sum(x)
  v<- sum((alpha.vec[1:2]-1)*log(x)) + (alpha.vec[3]-1)*log(xd)
  v1<- gamma(sum(alpha.vec))/prod(gamma(alpha.vec))*exp(v)
  v1  
}


#2 glmnet design matrix and reponse matrix
gen.designR <- function(p)
{
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
  return(list(x = x, y = y))  
}


#3 Construct Qmat 
gen.Qmat <- function(p,counter, tn, rand.sim = F )
{
  if (counter == 1)
  {
    alpha.vec1 = c(10, 8, 3)
    alpha.vec2 = c(3, 3, 10)
    alpha.vec3 = c(10, 3, 4)
  }
  
  if (counter == 2) 
  {
    alpha.vec1 = c(1, 8, 3)
    alpha.vec2 = c(3, 3, 10)
    alpha.vec3 = c(10, 3, 4)}
  
  if (counter == 3)
  {
    alpha.vec1 = c(3, 2, 1)
    alpha.vec2 = c(3, 4, 1)
    alpha.vec3 = c(3, 6, 1)
  }
  
  if (counter<4) 
  {
    pi.var1 = 1/3
    pi.var2 = 1/3
    pi.var3 = 1/3
    z = seq(1/p,1-1/p, length = p)
    dmat <- matrix(0, length(z), length(z))
    for (i in 1:length(z))
    {
      for (j in 1:length(z))
      {
        xvar = c(z[i], z[j])
        if (sum(xvar)>= 1) next
        dmat[i,j] <- pi.var1*ddirichlet(xvar,alpha.vec1)+
          pi.var2*ddirichlet(xvar,alpha.vec2) +
          (1-pi.var1-pi.var2)*ddirichlet(xvar,alpha.vec3) 
      }
    }
    #dmat
    dmat <- dmat/sum(dmat)
    if (rand.sim == T)
    {
      tmp <- drop(rmultinom(1, tn, prob = c(dmat)) )
      dmat1 <- matrix(tmp, p, p)
    }
    
    if (rand.sim == F) dmat1 <- round(dmat*tn)
    
    dmat2 <- matrix(0,p,p)
    for (j in 1:p) dmat2[,j] <- rev(dmat1[j,])
    dmat2 <- dmat2 + t(dmat2) 
    Qmat = dmat2
  }
  
  if (counter == 4)
  {
    #Qmat = matrix(  rpois(  p^2 , p^2)^2 , p, p)
    Qmat =  matrix(  runif(p^2), p , p)
    Qmat[lower.tri(Qmat, diag = FALSE)] = 0
    #Qmat[sample(1:length(Qmat), trunc(length(Qmat)/4))] = 0
    Qmat =  Qmat + t(Qmat)
    diag(Qmat) <- 0
    dmat1 = NULL
    n = sum(Qmat)
    nj = colSums(Qmat)
    Qpmat = Qmat
    for (j in 1:nrow(Qmat)) Qpmat[j,] = round(Qmat[j,]/n*2,2)
    return(list(Qmat = Qmat, dmat = dmat1, Qpmat = Qpmat))
    
  }
  
  n = sum(Qmat)
  nj = colSums(Qmat)
  Qpmat = Qmat
  for (j in 1:nrow(Qmat)) Qpmat[j,] = Qmat[j,]/n*2
  
    # dmat1 has the coordiantes for image()
    return(list(Qmat = Qmat, dmat = dmat1, Qpmat = Qpmat))
}
  


#4 simulation code in Gmat.hat
gen.Gmathat <- function(Gmat, Qmat)
{
  p = ncol(Gmat)
  gmat.prob<-c(Gmat)
  gmat.num <- c(Qmat)
  gmat.gen<- rep(0, length(gmat.num))
  for (i in 1:length(gmat.num))
  {
    gmat.gen[i] <- rbinom(n = 1, size = gmat.num[i], prob = gmat.prob[i])
  }
  
  Gmat.hat <- matrix(gmat.gen,p,p)
  Gmat.hat[lower.tri(Gmat.hat, diag = T)] = 0
  tmp <- Qmat - t(Gmat.hat)
  Gmat.hat[lower.tri(Qmat)]<- tmp[lower.tri(Qmat)]
  Gmat.hat <- Gmat.hat/Qmat
  Gmat.hat[!is.finite(Gmat.hat)] = 0
  return( Gmat.hat )
}

#5 naive est
naive.fun <- function(Qpmat, Gmat.hat,x,y,p)
{
  wmat = Qpmat*Gmat.hat
  wmat = t(wmat)
  wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))] 
  # fit glmnet
  fit <- glmnet(x, y, family = 'binomial',
                intercept = FALSE, weights = wvec, lambda = 0, standardize = F, thresh = 1e-09)
  est = c(fit$beta[,1],0)
  naive.est <- est
  return(naive.est)
}
##

# # check the existence of estimator by Qpmat.c1
# # note that the condition of existence only holds for population version
# # 1. make edges  
# tmp <- cbind( rep(1:p, p), rep(1:p, each = p), c(Qpmat.c1))
# tmp = tmp[tmp[,3] > 0, 1:2]
# ed = c(t(tmp))
# # 2. make directed graph
# g <-make_graph(ed, n = p, directed = TRUE)
# # 3. connectivity check
# cg <- components(g, mode = 'strong')
# # 4. report
# if (cg$no  ==  1) print("estimator exists!") else print("estimator does not exist !")
# 


prop.fun <- function(max.k, Qpmat, Gmat.hat, x, y, p)
{
  Result = NULL
  Result.list = list()
  k = 0
  for (k in 0:max.k)
  {
    #cat('outer interataion::: ', k , '\n')
    # Thresholding Qpmat.c1 by constant 'cvar'
    # unique cvec
    cvec<- sort(unique(Qpmat[upper.tri(Qpmat)]))
    # the first threshold
    if (k == 0 ) idx <- which(Qpmat <= 0)  
    if ( k>0 ) idx <- which(Qpmat <= cvec[k]) 
    Qpmat.c1 = Qpmat  
    Qpmat.c1[idx] <- 0
    ###############################
    # set weight-vector 
    i1 = 1 ; i2 = 2
    idx = 1
    result = matrix(0, p*(p-1)/2, 4)
    for ( i1 in 1:(p-1))
    {
      for (i2 in (i1+1):p) 
      {
        Qpmat.c2 = Qpmat.c1
        nvec1 = Qpmat.c1[i1,]
        nvec2 = Qpmat.c1[i2,]
        idx1 = which(nvec1 == 0 | nvec2 == 0)
        idx2  =  setdiff( idx1, c(i1, i2))
        nvec3 = (nvec1[-idx1]+nvec2[-idx1])/2
        Qpmat.c2[i1,-idx1] = Qpmat.c2[i2,-idx1] = nvec3
        if (length(idx2)>0)  Qpmat.c2[i1,idx2] =  Qpmat.c2[i2,idx2] = 0
        Qpmat.c2[,i1] <- Qpmat.c2[i1,]
        Qpmat.c2[,i2] <- Qpmat.c2[i2,]
        # check the result:  Qpmat.c2[i1,]
        #                 :  Qpmat.c2[i2,]
        
        # if (length(idx2) > 0 & max( c( Qpmat.c1[i1,idx2],Qpmat.c1[i2,idx2]) ) >0) :: check the existence of estimator 
        
        ###############################
        # set weight-vector in glmnet
        wmat = Qpmat.c2*Gmat.hat
        wmat = t(wmat)
        wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))] 
        # fit glmnet
        fit <- glmnet(x, y, family = 'binomial',
                      intercept = FALSE, weights = wvec, lambda = 0, standardize = F, thresh = 1e-09)
        
        est = c(fit$beta[,1],0)
        result[idx, 1:2] = c(i1, i2)
        if( est[i1] > est[i2]) result[idx, 3] = 1
        
        
        pmat = plogis(outer(est, est, FUN = '-'))
        vmat = pmat*(1-pmat)
        
        v1 = rowSums(Qpmat.c2*vmat)
        v1 <- v1[-p]
        v2 = (-Qpmat.c2*vmat)[-p,-p]
        diag(v2) <- v1
        inv.v2 <- solve(v2)
        if (i2 < p )  result[idx, 4] <- inv.v2[i1,i1] + inv.v2[i2,i2] - 2*inv.v2[i1,i2]
        if (i2 == p ) result[idx, 4] <- inv.v2[i1,i1]
        idx = idx + 1
        #cat('     inner iteration:' , idx, '\n')
        #plot(c(fit$beta[,1],0), type ='b')
        #min(-diff(est, 1))  
        #sum(-diff(est, 1)<0)
        #v1 = rowSums(Qpmat.c2*vmat)
      }
    }
    Result <- cbind(Result, result[,4])
    Result.list[[k+1]] <- result
  }
  return(list(Result = Result, Result.list = Result.list))
  
}


#6 aggregation
agg.fun <- function(prop.result)    
{
  min.vec<-apply(prop.result$Result, 1, which.min)
  agg.result <- matrix(0,nrow(prop.result$Result), 3)
  agg.result[,1:2] <- (prop.result$Result.list[[1]])[,1:2]
  for (i in 1:nrow(agg.result))  
  {
    agg.result[i,3] <- (prop.result$Result.list[[min.vec[i]]])[i,3]
  }
  return(agg.result)  
}


summary.fun<- function(tmp) 
{
  sum(tmp[,3])== nrow(tmp)
  sum(tmp[,3]==0)
  #tmp[tmp[,3] == 1,1:2]
  tmp1 <- tmp[,1:3]
  tmp2 <- tmp[,c(2,1,3)]
  tmp2[,3]<- abs(tmp2[,3]-1)
  sc <- c()
  for ( i in 1:p)
  {
    sc[i] = sum(tmp1[tmp1[,1]== i ,3]) + sum(tmp2[tmp2[,1]== i ,3])
  }
  return(sc)  
}



###
cons.rank<- function(tmp)
{
  
  tmp.copy<- tmp[,c(2,1,3)]
  tmp.copy[,3] <- 1-tmp.copy[,3]
  tmp <- rbind(tmp, tmp.copy)
  tmp <- tmp[tmp[,1] != 0,]
  
  Cset <- unique( tmp[,1])
  sc <- c()
  while( length(Cset)>0)
  {
    Uset <- c()
    Cset <- unique( tmp[,1])
    for (i in Cset)
    {
      tmp.sub1<-tmp[tmp[,1] == i,, drop = F]
      if ( sum(tmp.sub1[,3]) == nrow(tmp.sub1) ) Uset <- c(Uset, i)
    }
    
    #cat(sc,'\n')
    if ( length(Uset) == 1)  
    {
      sc <- c(sc,Uset)
      if (length(sc)==(p-1)) 
      {
        sc<- c( sc, setdiff(Cset, Uset) )
        break
      }
      Uset.idx<-tmp[, 1] == Uset | tmp[, 2] == Uset
      tmp <- tmp[!Uset.idx,]
    } else {
      sc <- NA
      break
    }
    
  }
  return(sc)
}


## cv code 

# make cv_mat
# 
cv_mat_fun <- function(Gmat_obs, Qmat, k_fold) 
{
  cv_mat = matrix(0, tn, 4)  
  s1 = 1
  for (j in 1:(ncol(Gmat_obs)-1))
  {
    for(i in (j+1):nrow(Gmat_obs))
    {
      a1 = Qmat[i,j]
      if (a1 == 0 ) next
      f_idx <- s1:(s1+a1-1)
      cv_mat[f_idx,1] <- i
      cv_mat[f_idx,2] <- j
      s1 = s1 + a1
      a2 = Gmat_obs[i,j]
      if (a2 == 0 ) next
      f_idx2 <- sample(f_idx,a2)
      cv_mat[f_idx2,3] = 1
    }
  }
  cv_mat[,4] <- sample(1:k_fold, tn, replace = TRUE)
  colnames(cv_mat) = c("j", "k", "y_jk", "partition")
  return(cv_mat)
}

cv_table_fun = function(cv_m)
{
  cv_m <- as.data.frame(cv_m)
  # require(dplyr)
  result = cv_m %>% group_by(j, k) %>% 
    summarize(sum = length(y_jk)) %>% as.matrix()
  
  Qmat_tr = matrix(0, p, p)
  for ( i in 1:nrow(result))
  {
    Qmat_tr[result[i, 1], result[i, 2]] <- result[i, 3]  
  }
  Qmat_tr <- Qmat_tr + t(Qmat_tr)
  
  result = cv_m %>% group_by(j, k) %>% 
    summarize(sum = sum(y_jk)) %>% as.matrix()
  
  Gmat_tr = matrix(0, p, p)
  for ( i in 1:nrow(result))
  {
    Gmat_tr[result[i, 1], result[i, 2]] <- result[i, 3]        
  }
  tmp <- Qmat_tr - t(Gmat_tr)
  Gmat_tr[upper.tri(Qmat_tr)] <- tmp[upper.tri(Qmat_tr)]
  return(list(Q=Qmat_tr, G=Gmat_tr))
}


gen_sim_fun = function(Gmat, Qmat)
{
  ## Gmat.hat : true Gmat을 이용하여 data generation을 할 경우에 승패 수와 전체 대결 수를 이용하여 
  ##            만든 Gmat의 추정값 
  gmat_prob <- c(Gmat)  ## Gmat_jk : j object와 k object에서 j가 k를 이길 확률.
  gmat_num <- c(Qmat)
  gmat_gen<- rep(0, length(gmat_num))
  for (i in 1:length(gmat_num))
  {
    gmat_gen[i] <- rbinom(n = 1, size = gmat_num[i], prob = gmat_prob[i])
  }
  Gmat_obs <- matrix(gmat_gen,p,p)
  Gmat_obs[lower.tri(Gmat_obs, diag = T)] = 0
  tmp <- Qmat - t(Gmat_obs)
  Gmat_obs[lower.tri(Qmat)]<- tmp[lower.tri(Qmat)]
  return( list(G = Gmat_obs, Q = Qmat) )
}


bt_fun = function(gen_fit, lambda.vec = NULL)
{ 
  Gmat.hat <- gen_fit$G
  Qmat <- gen_fit$Q
  p = ncol(Qmat)
  Gmat.hat <- Gmat.hat/Qmat
  Gmat.hat[!is.finite(Gmat.hat)] = 0
  n = sum(Qmat)
  Qpmat = Qmat/n*2
  wmat = Qpmat*Gmat.hat  
  wmat = t(wmat)
  wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  
  # fit glmnet
  fit <- glmnet(x, y, family = 'binomial',
                intercept = FALSE, weights = wvec, lambda = 0, 
                standardize = F, thresh = 1e-09)
  est = c(fit$beta[,1],0)
  if (is.null(lambda.vec)) cor.r = NULL else cor.r = cor(est, lambda.vec, method = 'kendall')    
  return( list (coefficients = est,
                cor = cor.r) )
  
}



gbt_step1_fun = function(Qpmat, Gmat.hat, p, cval)      
{
  result = matrix(0, p*(p-1)/2, 4)
  idx = 1
  # select a pair of items for obtaining rank consistent estimator
  # define a matrix to save paired results
  for ( i1 in 1:(p-1))
  {
    for (i2 in (i1+1):p)
    {
      Qpmat.c1 = Qpmat
      # threshold step
      idx1 <- ( Qpmat.c1[i1,] <= cval )      
      idx2 <- ( Qpmat.c1[i2,] <= cval )  ## intersect(!idx1,!idx2)=\hat{O}_jk   
      if (sum(idx1)>0 )  ##length->sum으로 고침 
      {
        Qpmat.c1[i1,idx1] <- 0 ;  Qpmat.c1[idx1,i1] <- 0
      }
      
      if (sum(idx2)>0 )  ##length->sum으로 고침
      {
        Qpmat.c1[i2,idx2] <- 0 ;  Qpmat.c1[idx2,i2] <- 0
      }
      
      Qpmat.c2 = Qpmat.c1
      ## thresholding procedure
      Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cval)  
      
      nvec1 = Qpmat.c1[i1,]
      nvec2 = Qpmat.c1[i2,]
      idx1 = which(nvec1 == 0 | nvec2 == 0) ## !idx1 : \hat{O}_jk
      idx2  =  setdiff( idx1, c(i1, i2))    
      
      # balancing the weight parameter in gBT model for the selected pair
      nvec3 = (nvec1[-idx1]+nvec2[-idx1])/2
      Qpmat.c2[i1,-idx1] = Qpmat.c2[i2,-idx1] = nvec3
      
      if (length(idx2)>0)  Qpmat.c2[i1,idx2] =  Qpmat.c2[i2,idx2] = 0
      
      Qpmat.c2[,i1] <- Qpmat.c2[i1,]  ## 대칭이 되도록 만들자 
      Qpmat.c2[,i2] <- Qpmat.c2[i2,]  ## 대칭이 되도록 만들자
      
      ## find V_jk(maximum connected set)                
      i1i2_adj_matrix = matrix(as.integer(Qpmat.c2>0) , p , p)  ## adjacency matrix
      i1i2_graph = graph_from_adjacency_matrix(i1i2_adj_matrix , mode="undirected" , weighted=NULL) ## make a graph
      i1i2_clusters = clusters(i1i2_graph)$mem ## clustering using adj matrix
      if (i1i2_clusters[i1] != i1i2_clusters[i2])
      {
        ## i1과 i2가 다른 connected 되지 않은 경우
        #cat('   k:',k,', ',i1,'and',i2, 'is not connected!!\n') 
        idx = idx + 1
        next  
      } 
      ## idx3 : edge index set of V_jk
      idx3 = sort(which(i1i2_clusters %in% i1i2_clusters[i1])) 
      
      #########################################
      ## computing gBT estimator
      #########################################
      wmat <- Qpmat.c2[idx3, idx3]*Gmat.hat[idx3, idx3]
      wmat = t(wmat)
      pp <- length(idx3)
      wvec = wmat[ - (1 + ( 0:(pp-1) ) *(pp+1))]  ## w_jj는 제거.. 
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
      
      # note that the gBT estimator may not exist because of the thresolding
      # use ridge
      try.fit <- try(fit <- glmnet(xx, yy, family = 'binomial',
                                   intercept = FALSE, weights = wvec, 
                                   lambda = 1e-5, alpha = 0, standardize = F, 
                                   thresh = 1e-09), silent = T)
      
      if (class(try.fit)[1] == 'try-error')
      {
        ## Q: result[idx,]에 error가 났다는 어떤 표시도 하지 않나..?? 
        ## A: result에 (i1,i2) 대신 (0,0)을 표시.. 
        idx = idx + 1
        next
      }
      
      
      est = c(fit$beta[,1],0) ## lambda_pp 추가 
      result[idx, 1:2] = c(i1, i2)
      # compare the values of lambda_{i1} and lambda_{i2}
      if( est[which(idx3==i1)] > est[which(idx3==i2)]) result[idx, 3] = 1
      ## \hat{lambda_i1}>\hat{lambda_i2}면 result[idx,3]에 1을 부여..(아니면 0) 
      
      
      ## calculate weight v_jk
      
      # 1. obtain the asymptotic variance
      pmat = plogis(outer(est, est, FUN = '-'))
      vmat = pmat*(1-pmat)
      v1 = rowSums(Qpmat.c2[idx3,idx3]*vmat)
      v1 <- v1[-pp]
      v2 = (-Qpmat.c2[idx3,idx3]*vmat)[-pp,-pp]
      diag(v2) <- v1
      ## calculate asymptotic covariance matrix
      tri.fit2 = try(inv.v2 <- solve(v2*sum(Qmat)/2))
      if (class(tri.fit2)[1] == 'try-error'){
        cat('   k:',k,', ',i1,'and',i2, ': cannot calaulate inverse')
        result[idx , 4] = 4/(sum(Qpmat.c2[idx3,idx3]*sum(Qmat)/2)/2)  ## alternative v_jk 
        idx = idx+1 ## if error -> next step
        next
      }
      i1_ind = which(idx3==i1); i2_ind = which(idx3==i2)
      if ((i1_ind<pp) & (i2_ind<pp))
      {
        result[idx, 4] <- inv.v2[i1_ind,i1_ind] + inv.v2[i2_ind,i2_ind] - 
          2*inv.v2[i1_ind,i2_ind]
      }
      
      if ((i1_ind==pp) | (i2_ind==pp) )
      {
        min_ind = min(i1_ind , i2_ind)
        result[idx , 4] = inv.v2[min_ind , min_ind]
      } 
      
      ## 2. v_jk = \sum n_ml
      #result[idx , 4] = sum(Qpmat.c2[idx3,idx3])/2
      
      #cat('     inner iteration:' , idx, '\n')
      idx = idx + 1
      
      #plot(c(fit$beta[,1],0), type ='b')
      #min(-diff(est, 1))  
      #sum(-diff(est, 1)<0)
      #v1 = rowSums(Qpmat.c2*vmat)
      #}
      #}
    }
  }
  return(result)  
}

# note that gbt_step2_fun has two types of returns: 
gbt_step2_fun = function(result, p, lambda.vec, newdata = NULL)
{
  tmp<-result
  not0_ind = (tmp[,1]!=0)
  tmp <-tmp[not0_ind, 1:3]
  p.set <-sort(unique(c(tmp[,1:2])))
  
  if (length(p.set) != p) 
  {
    return(list(gbt_est = rep(NA,p),
                cor = NA, test_cor = NA))
  }
  
  xx <- matrix(0, nrow(tmp)*2, p)
  yy <- rep(0, nrow(tmp)*2)
  i = 1
  for ( i in 1:nrow(tmp))
  {
    vec1<-tmp[i,1:2]; vec2<- tmp[i,3]
    xx[2*(i-1)+1, vec1] <- c(1,-1) ; yy[2*(i-1)+1] <- vec2
    xx[2*i, vec1] <- c(-1,1) ; yy[2*i] <- abs(vec2 - 1)
  }
  xx<- xx[,-p]
  
  ## See fit: Note that weights denotes v_jk 
  fit<-glmnet(xx,yy, family = 'binomial', alpha = 0, lambda = 1e-5, 
              intercept = FALSE,
              weights = rep(result[not0_ind, 4],each=2) , standardize = F)
  
  gbt.est <- c(fit$beta[,1],0)
  
  if (is.null(lambda.vec)) 
  {
    cor.r = NA
  } else {
    cor.r <- cor(gbt.est, lambda.vec, method = 'kendall')    
  }
  
  if (is.null(newdata))
  {
    test_cor = NA
  } else {
    tmp = matrix(0,p,p)
    Q = newdata$Q
    G = newdata$G
    for (i in 1:p)
    {
      for (j in 1:p)
      {
        if (i == j) next
        tmp[i,j] = as.integer( gbt.est[i]-gbt.est[j] > 0 )
      }
    }
    test_cor = sum(tmp*G)/(sum(Q)/2)
  }
  return(list(gbt_est = gbt.est,
      cor = cor.r, test_cor = test_cor))
}



sparse_gen_fun <- function(dmat, kn, rn, tn)
{
  dmat1 <- dmat ## dmat : {q_jk} matrix
  u.idx <- which( dmat > 0) ## q_jk가 0보다 큰 index 
  sel.u.idx<- sample(u.idx, kn) ## d개를 sampling함. 
  dmat1[sel.u.idx]  <- 0  ## d개의 선택된 q_jk에는 0을 대입(어차피 해당 n_jk은 n_s로 고정할 것이므로) 
  dmat1 <- dmat1/sum(dmat1) ## dmat1을 prob distribution으로 만들어주기 위해 normalize. 
  d.sample <- drop (rmultinom(1, tn-rn*kn, prob = c(dmat1)))  ## n_jk 만들기 
  d.sample[sel.u.idx] <- rn ## d개의 선택된 n_jk에 n_s를 대입 
  dmat1 <- matrix(d.sample, p , p)  ## matrix 형태로 만들어줌. 
  Qmat <- matrix(0, p, p )
  for (j in 1:p) Qmat[,j] <- rev(dmat1[j,])
  Qmat <- Qmat + t(Qmat)  ## Qmat : 최종적인 n_jk matrix 
  return(Qmat)
}


cv.gbt_fun  = function(gen_fit, cvec, k_fold, lambda.vec)
{
  fit = gen_fit
  p = ncol(fit$Q)
  cv_mat = cv_mat_fun(fit$G, fit$Q, k_fold) 
  cor_mat = matrix(NA, k_fold, length(cvec))
  for (k in 1:length(cvec))
  {
    ######### gBT model ###########
    cval <- cvec[k]
    # k-fold
    for (k_num in 1:k_fold)
    {
      result = NULL
      tmp_te = cv_mat[cv_mat[,"partition"] == k_num,-4]
      tmp_tr = cv_mat[cv_mat[,"partition"] != k_num,-4]
      cv_table <- cv_table_fun(tmp_tr)
      Qmat <- cv_table$Q
      Gmat.hat <- cv_table$G
      Gmat.hat = Gmat.hat/Qmat
      Gmat.hat[!is.finite(Gmat.hat)] = 0
      n = sum(Qmat)
      Qpmat = Qmat/n*2
      result <- gbt_step1_fun(Qpmat, Gmat.hat, p, cval)
      # gbt_step2_fun
      cv_table <- cv_table_fun(tmp_te)
      gbt_fit <- gbt_step2_fun(result, p, lambda.vec, cv_table)
#      cat ("k:", k, "  k_num:", k_num,"  ", gbt_fit$test_cor,'\n')
      if (any(is.na(gbt_fit))) next 
      cor_mat[k_num, k] <- gbt_fit$test_cor
    }
  }
  return(cor_mat)
}


gbt_fun = function(gen_fit, cval, lambda.vec)
{
  Gmat.hat <- gen_fit$G
  Qmat <- gen_fit$Q
  p = ncol(Qmat)
  Gmat.hat <- Gmat.hat/Qmat
  Gmat.hat[!is.finite(Gmat.hat)] = 0
  n = sum(Qmat)
  Qpmat = Qmat/n*2
  result <- gbt_step1_fun(Qpmat, Gmat.hat, p, cval)
  gbt_fit<- gbt_step2_fun(result, p, lambda.vec) 
  gbt_est = gbt_fit$gbt_est
  cor = gbt_fit$cor
  return( list (coefficients = gbt_est,
                cor = cor) )
}