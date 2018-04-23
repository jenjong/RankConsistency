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
