## description
## investigate the variance of the proposed estimator
rm(list = ls()); gc()
#setwd("E:\\rank consistency\\simulation\\code\\")
#setwd("C:/Users/uos_stat/Dropbox/A rank-consistency/prog/temp")
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
require('glmnet')
library(igraph)
library(ggplot2)
source('./lib/sim.R')

##########################################################################################
max.k = 10
p = 10
kn <- 10  ## d
rn <- 1   ## n_s
df = 1    ## t분포 자유도 
counter = 1
sim.iter = 200    ## 전체 simulation 과정 반복 수 
source('./lib/exe-2.R')
tn_vec = c(500,5000,50000)
cor.naive_list = list()
cor.r_list = list()
k_fold = 5
tn_i = 1
for (tn_i in 1:3)
{
  tn = tn_vec[tn_i]  ## tn 정의 (전체 rank pair의 수.)
  cat ('total n:' , tn , '\n')
  cor.naive<- rep(0,sim.iter) ## BT를 이용한 kendall's tau 저장하는 벡터 
  cor.r <- matrix(0,sim.iter,max.k) ## gBT를 이용한 kendall's tau 저장하는 벡터 
  ii = 1
  for  (ii in 1:sim.iter)
  {
    if (ii %% 10 == 0)  cat(' ',ii,'-th iteration\n')
    
    set.seed(ii+123) ## set seed
    
    ### generating number of comparison using dmat
    ### output: Qmat (n_jk matrix)
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
    
    cvec<-(0:(max.k-1))/tn  ## cvec : threshold c들이 모여있는 벡터 
    
    ##############################
    # generating the result of win-loss ratio using Gmat(true probability matrix) on Qmat (number of matching)
    # output: Gmat.hat
    ## Gmat.hat : true Gmat을 이용하여 data generation을 할 경우에 승패 수와 전체 대결 수를 이용하여 
    ##            만든 Gmat의 추정값 
    gmat.prob<-c(Gmat)  ## Gmat_jk : j object와 k object에서 j가 k를 이길 확률.
    ## 논문에서는 F*를 t분포라 했는데 여기서는 왜
    ## G를 t분포라 한건지..?? 
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
    cv_mat = matrix(0, tn, 4)
    
    s1 = 1
    for (j in 1:(ncol(Gmat)-1))
    {
      for(i in (j+1):nrow(Gmat))
      {
        a1 = Qmat[i,j]
        if (a1 == 0 ) next
        f_idx <- s1:(s1+a1-1)
        cv_mat[f_idx,1] <- i
        cv_mat[f_idx,2] <- j
        s1 = s1 + a1
        a2 = Gmat.hat[i,j]
        if (a2 == 0 ) next
        f_idx2 <- sample(f_idx,a2)
        cv_mat[f_idx,3] = 1
      }
    }
    cv_mat[,4] <- sample(1:k_fold, tn, replace = TRUE)
    cv_te = cv_mat[cv_mat[,4] == 1,]
    cv_tr = cv_mat[cv_mat[,4] != 1,]
    Qmat_tr = matrix(0, )
    Gmat_tr =
    for ()
    
    
    
    Gmat.hat <- Gmat.hat/Qmat
    Gmat.hat[!is.finite(Gmat.hat)] = 0
    
    ###############################
    # Generating Qpmat that consists of q_{jk}
    n = sum(Qmat)     ## n : 2*tn (total n의 2배)
    Qpmat = Qmat/n*2  ## 왜 2를 곱해주지...??? -> n이 tn의 2배이므로..(Qpmat=Qmat/tn) 
    
    ###############################
    # define the variable to restore the simulation results
    Result = NULL
    Naive.list = list()
    Result.list = list()
    #cvec<- sort(unique(Qpmat[upper.tri(Qpmat)]))
    
    ################################
    # Start the algrotihm for each c (threshold variable)
    for (k in 0:(max.k-1))
    {
      #cat ('thresholding value c:',cvec[k+1],'\n')
      #cat('outer interataion::: ', k , '\n')
      Qpmat.c1 = Qpmat
      
      ##############################
      ######## BT model ############
      # output: naive.est
      if (k == 0)
      {
        wmat = Qpmat.c1*Gmat.hat  ## j가 k를 이긴 횟수에 비례해서 weight를 만든다. 
        wmat = t(wmat)
        wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  ## w_jj는 제거.. 
        ## 일반적인 BT model을 fit하자.. 
        # fit glmnet
        fit <- glmnet(x, y, family = 'binomial',
                      intercept = FALSE, weights = wvec, lambda = 0, standardize = F, thresh = 1e-09)
        est = c(fit$beta[,1],0) ## lambda hat들을 저장. (lambda_p = 0)
        naive.est <- est        ## naive.est에 est를 저장. 
      }
      
      ###############################
      ######### gBT model ###########
      # set weight-vector
      # select a pair of items for obtaining rank consistent estimator
      idx = 1
      # define a matrix to save paired results
      result = matrix(0, p*(p-1)/2, 4)
      for ( i1 in 1:(p-1))
      {
        for (i2 in (i1+1):p)
        {
          Qpmat.c1 = Qpmat
          # threshold step
          idx1 <- ( Qpmat.c1[i1,] <= cvec[k+1] )      
          idx2 <- ( Qpmat.c1[i2,] <= cvec[k+1] )  ## intersect(!idx1,!idx2)=\hat{O}_jk   
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
          Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cvec[k+1])  
          
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
          if (i1i2_clusters[i1] != i1i2_clusters[i2]){  ## i1과 i2가 다른 connected 되지 않은 경우
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
                                       intercept = FALSE, weights = wvec, lambda = 1e-5, alpha = 0, standardize = F, thresh = 1e-09),
                         silent = T)
          
          if (class(try.fit)[1] == 'try-error')
          {
            ## Q: result[idx,]에 error가 났다는 어떤 표시도 하지 않나..?? 
            ## A: result에 (i1,i2) 대신 (0,0)을 표시.. 
            idx = idx + 1
            next
          }
          
          # compare the values of lambda_{i1} and lambda_{i2}
          est = c(fit$beta[,1],0) ## lambda_pp 추가 
          result[idx, 1:2] = c(i1, i2)
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
          if ((i1_ind<pp) & (i2_ind<pp)){
            result[idx, 4] <- inv.v2[i1_ind,i1_ind] + inv.v2[i2_ind,i2_ind] - 
              2*inv.v2[i1_ind,i2_ind]
          }
          if ((i1_ind==pp) | (i2_ind==pp) ){
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
      
      Result <- cbind(Result, result[,4]) ## Result는 뭐지..?? 
      ## result[,4]는 0만 있는 vector들임..저장하는 이유가? 
      Result.list[[k+1]] <- result
    }
    
    #############################################################################
    
    cor.naive[ii] <- cor(naive.est, lambda.vec, method = 'kendall') 
    ## BTmodel을 이용한 결과에 대한 kendall's tau 저장. 
    
    for (k in 1:max.k)
    {
      tmp<-Result.list[[k]]
      not0_ind = (tmp[,1]!=0)
      tmp <-tmp[not0_ind, 1:3]  ## glm 잘 돌아간 애들만 뽑는 과정..
      p.set <-sort(unique(c(tmp[,1:2])))
      
      # corrrrrrr~!!
      if (length(p.set) != p) ## not rank-recoverable..?? 맞나..?? 
      {
        cor.r[ii,k] <- NA
        next
      }
      
      ##### gBT refitting
      
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
      
      fit<-glmnet(xx,yy, family = 'binomial', alpha = 0, lambda = 1e-5, intercept = FALSE,
                  weights = rep(Result[not0_ind,k],each=2) , standardize = F)
      ## weight vector v_jk는 들어가지 않나?? 
      gbt.est <- c(fit$beta[,1],0)
      cor.r[ii,k] <- cor(gbt.est, lambda.vec, method = 'kendall')
      
    }
    # end of max.k iter
    #cat( which.max( colMeans(cor.r[1:ii, , drop = F], na.rm = T)), '\n')
    #cat(which.max( colMeans(cor.topo[1:ii,, drop = F], na.rm = T)), '\n')
    #cat( mean(cor.naive[1:ii]), '\n')    
  }
  cor.naive_list[[tn_i]] = cor.naive
  cor.r_list[[tn_i]] = cor.r
}

cor.r_list1 = cor.r_list
cor.naive_list1 = cor.naive_list


##########################################################################################
rn <- 3   ## n_s

cor.naive_list = list()
cor.r_list = list()

for (tn_i in 1:3){
  tn = tn_vec[tn_i]  ## tn 정의 
  cat ('total n:' , tn , '\n')
  cor.naive<- rep(0,sim.iter) ## BT를 이용한 kendall's tau 저장하는 벡터 
  cor.r <- matrix(0,sim.iter,max.k) ## gBT를 이용한 kendall's tau 저장하는 벡터 
  
  for  (ii in 1:sim.iter)
  {
    if (ii %% 10 == 0){
      cat(' ',ii,'-th iteration\n')
    }
    set.seed(ii+123) ## set seed
    
    ### generating number of comparison using dmat
    ### output: Qmat (n_jk matrix)
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
    
    cvec<-(0:(max.k-1))/tn  ## cvec : threshold c들이 모여있는 벡터 
    
    ##############################
    # generating the result of win-loss ratio using Gmat(true probability matrix) on Qmat (number of matching)
    # output: Gmat.hat
    ## Gmat.hat : true Gmat을 이용하여 data generation을 할 경우에 승패 수와 전체 대결 수를 이용하여 
    ##            만든 Gmat의 추정값 
    gmat.prob<-c(Gmat)  ## Gmat_jk : j object와 k object에서 j가 k를 이길 확률.
    ## 논문에서는 F*를 t분포라 했는데 여기서는 왜
    ## G를 t분포라 한건지..?? 
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
    
    ###############################
    # Generating Qpmat that consists of q_{jk}
    n = sum(Qmat)     ## n : 2*tn (total n의 2배)
    Qpmat = Qmat/n*2  ## 왜 2를 곱해주지...??? -> n이 tn의 2배이므로..(Qpmat=Qmat/tn) 
    
    ###############################
    # define the variable to restore the simulation results
    Result = NULL
    Naive.list = list()
    Result.list = list()
    #cvec<- sort(unique(Qpmat[upper.tri(Qpmat)]))
    
    ################################
    # Start the algrotihm for each c (threshold variable)
    for (k in 0:(max.k-1))
    {
      #cat ('thresholding value c:',cvec[k+1],'\n')
      #cat('outer interataion::: ', k , '\n')
      Qpmat.c1 = Qpmat
      
      ##############################
      ######## BT model ############
      # output: naive.est
      if (k == 0)
      {
        wmat = Qpmat.c1*Gmat.hat  ## j가 k를 이긴 횟수에 비례해서 weight를 만든다. 
        wmat = t(wmat)
        wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  ## w_jj는 제거.. 
        ## 일반적인 BT model을 fit하자.. 
        # fit glmnet
        fit <- glmnet(x, y, family = 'binomial',
                      intercept = FALSE, weights = wvec, lambda = 0, standardize = F, thresh = 1e-09)
        est = c(fit$beta[,1],0) ## lambda hat들을 저장. (lambda_p = 0)
        naive.est <- est        ## naive.est에 est를 저장. 
      }
      
      ###############################
      ######### gBT model ###########
      # set weight-vector
      # select a pair of items for obtaining rank consistent estimator
      idx = 1
      # define a matrix to save paired results
      result = matrix(0, p*(p-1)/2, 4)
      for ( i1 in 1:(p-1))
      {
        for (i2 in (i1+1):p)
        {
          Qpmat.c1 = Qpmat
          # threshold step
          idx1 <- ( Qpmat.c1[i1,] <= cvec[k+1] )      
          idx2 <- ( Qpmat.c1[i2,] <= cvec[k+1] )  ## intersect(!idx1,!idx2)=\hat{O}_jk   
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
          Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cvec[k+1])  
          
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
          if (i1i2_clusters[i1] != i1i2_clusters[i2]){  ## i1과 i2가 다른 connected 되지 않은 경우
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
                                       intercept = FALSE, weights = wvec, lambda = 1e-5, alpha = 0, standardize = F, thresh = 1e-09),
                         silent = T)
          
          if (class(try.fit)[1] == 'try-error')
          {
            ## Q: result[idx,]에 error가 났다는 어떤 표시도 하지 않나..?? 
            ## A: result에 (i1,i2) 대신 (0,0)을 표시.. 
            idx = idx + 1
            next
          }
          
          # compare the values of lambda_{i1} and lambda_{i2}
          est = c(fit$beta[,1],0) ## lambda_pp 추가 
          result[idx, 1:2] = c(i1, i2)
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
          if ((i1_ind<pp) & (i2_ind<pp)){
            result[idx, 4] <- inv.v2[i1_ind,i1_ind] + inv.v2[i2_ind,i2_ind] - 
              2*inv.v2[i1_ind,i2_ind]
          }
          if ((i1_ind==pp) | (i2_ind==pp) ){
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
      
      Result <- cbind(Result, result[,4]) ## Result는 뭐지..?? 
      ## result[,4]는 0만 있는 vector들임..저장하는 이유가? 
      Result.list[[k+1]] <- result
    }
    
    #############################################################################
    
    cor.naive[ii] <- cor(naive.est, lambda.vec, method = 'kendall') 
    ## BTmodel을 이용한 결과에 대한 kendall's tau 저장. 
    
    for (k in 1:max.k)
    {
      tmp<-Result.list[[k]]
      not0_ind = (tmp[,1]!=0)
      tmp <-tmp[not0_ind, 1:3]  ## glm 잘 돌아간 애들만 뽑는 과정..
      p.set <-sort(unique(c(tmp[,1:2])))
      
      # corrrrrrr~!!
      if (length(p.set) != p) ## not rank-recoverable..?? 맞나..?? 
      {
        cor.r[ii,k] <- NA
        next
      }
      
      ##### gBT refitting
      
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
      
      fit<-glmnet(xx,yy, family = 'binomial', alpha = 0, lambda = 1e-5, intercept = FALSE,
                  weights = rep(Result[not0_ind,k],each=2) , standardize = F)
      ## weight vector v_jk는 들어가지 않나?? 
      gbt.est <- c(fit$beta[,1],0)
      cor.r[ii,k] <- cor(gbt.est, lambda.vec, method = 'kendall')
      
    }
    # end of max.k iter
    #cat( which.max( colMeans(cor.r[1:ii, , drop = F], na.rm = T)), '\n')
    #cat(which.max( colMeans(cor.topo[1:ii,, drop = F], na.rm = T)), '\n')
    #cat( mean(cor.naive[1:ii]), '\n')    
  }
  cor.naive_list[[tn_i]] = cor.naive
  cor.r_list[[tn_i]] = cor.r
}

cor.r_list3 = cor.r_list
cor.naive_list3 = cor.naive_list


##########################################################################################
rn <- 5   ## n_s

cor.naive_list = list()
cor.r_list = list()

for (tn_i in 1:3){
  tn = tn_vec[tn_i]  ## tn 정의 
  cat ('total n:' , tn , '\n')
  cor.naive<- rep(0,sim.iter) ## BT를 이용한 kendall's tau 저장하는 벡터 
  cor.r <- matrix(0,sim.iter,max.k) ## gBT를 이용한 kendall's tau 저장하는 벡터 
  
  for  (ii in 1:sim.iter)
  {
    if (ii %% 10 == 0){
      cat(' ',ii,'-th iteration\n')
    }
    set.seed(ii+123) ## set seed
    
    ### generating number of comparison using dmat
    ### output: Qmat (n_jk matrix)
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
    
    cvec<-(0:(max.k-1))/tn  ## cvec : threshold c들이 모여있는 벡터 
    
    ##############################
    # generating the result of win-loss ratio using Gmat(true probability matrix) on Qmat (number of matching)
    # output: Gmat.hat
    ## Gmat.hat : true Gmat을 이용하여 data generation을 할 경우에 승패 수와 전체 대결 수를 이용하여 
    ##            만든 Gmat의 추정값 
    gmat.prob<-c(Gmat)  ## Gmat_jk : j object와 k object에서 j가 k를 이길 확률.
    ## 논문에서는 F*를 t분포라 했는데 여기서는 왜
    ## G를 t분포라 한건지..?? 
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
    
    ###############################
    # Generating Qpmat that consists of q_{jk}
    n = sum(Qmat)     ## n : 2*tn (total n의 2배)
    Qpmat = Qmat/n*2  ## 왜 2를 곱해주지...??? -> n이 tn의 2배이므로..(Qpmat=Qmat/tn) 
    
    ###############################
    # define the variable to restore the simulation results
    Result = NULL
    Naive.list = list()
    Result.list = list()
    #cvec<- sort(unique(Qpmat[upper.tri(Qpmat)]))
    
    ################################
    # Start the algrotihm for each c (threshold variable)
    for (k in 0:(max.k-1))
    {
      #cat ('thresholding value c:',cvec[k+1],'\n')
      #cat('outer interataion::: ', k , '\n')
      Qpmat.c1 = Qpmat
      
      ##############################
      ######## BT model ############
      # output: naive.est
      if (k == 0)
      {
        wmat = Qpmat.c1*Gmat.hat  ## j가 k를 이긴 횟수에 비례해서 weight를 만든다. 
        wmat = t(wmat)
        wvec = wmat[ - (1 + ( 0:(p-1) ) *(p+1))]  ## w_jj는 제거.. 
        ## 일반적인 BT model을 fit하자.. 
        # fit glmnet
        fit <- glmnet(x, y, family = 'binomial',
                      intercept = FALSE, weights = wvec, lambda = 0, standardize = F, thresh = 1e-09)
        est = c(fit$beta[,1],0) ## lambda hat들을 저장. (lambda_p = 0)
        naive.est <- est        ## naive.est에 est를 저장. 
      }
      
      ###############################
      ######### gBT model ###########
      # set weight-vector
      # select a pair of items for obtaining rank consistent estimator
      idx = 1
      # define a matrix to save paired results
      result = matrix(0, p*(p-1)/2, 4)
      for ( i1 in 1:(p-1))
      {
        for (i2 in (i1+1):p)
        {
          Qpmat.c1 = Qpmat
          # threshold step
          idx1 <- ( Qpmat.c1[i1,] <= cvec[k+1] )      
          idx2 <- ( Qpmat.c1[i2,] <= cvec[k+1] )  ## intersect(!idx1,!idx2)=\hat{O}_jk   
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
          Qpmat.c2 = Qpmat.c2*(Qpmat.c2>cvec[k+1])  
          
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
          if (i1i2_clusters[i1] != i1i2_clusters[i2]){  ## i1과 i2가 다른 connected 되지 않은 경우
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
                                       intercept = FALSE, weights = wvec, lambda = 1e-5, alpha = 0, standardize = F, thresh = 1e-09),
                         silent = T)
          
          if (class(try.fit)[1] == 'try-error')
          {
            ## Q: result[idx,]에 error가 났다는 어떤 표시도 하지 않나..?? 
            ## A: result에 (i1,i2) 대신 (0,0)을 표시.. 
            idx = idx + 1
            next
          }
          
          # compare the values of lambda_{i1} and lambda_{i2}
          est = c(fit$beta[,1],0) ## lambda_pp 추가 
          result[idx, 1:2] = c(i1, i2)
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
          if ((i1_ind<pp) & (i2_ind<pp)){
            result[idx, 4] <- inv.v2[i1_ind,i1_ind] + inv.v2[i2_ind,i2_ind] - 
              2*inv.v2[i1_ind,i2_ind]
          }
          if ((i1_ind==pp) | (i2_ind==pp) ){
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
      
      Result <- cbind(Result, result[,4]) ## Result는 뭐지..?? 
      ## result[,4]는 0만 있는 vector들임..저장하는 이유가? 
      Result.list[[k+1]] <- result
    }
    
    #############################################################################
    
    cor.naive[ii] <- cor(naive.est, lambda.vec, method = 'kendall') 
    ## BTmodel을 이용한 결과에 대한 kendall's tau 저장. 
    
    for (k in 1:max.k)
    {
      tmp<-Result.list[[k]]
      not0_ind = (tmp[,1]!=0)
      tmp <-tmp[not0_ind, 1:3]  ## glm 잘 돌아간 애들만 뽑는 과정..
      p.set <-sort(unique(c(tmp[,1:2])))
      
      # corrrrrrr~!!
      if (length(p.set) != p) ## not rank-recoverable..?? 맞나..?? 
      {
        cor.r[ii,k] <- NA
        next
      }
      
      ##### gBT refitting
      
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
      
      fit<-glmnet(xx,yy, family = 'binomial', alpha = 0, lambda = 1e-5, intercept = FALSE,
                  weights = rep(Result[not0_ind,k],each=2) , standardize = F)
      ## weight vector v_jk는 들어가지 않나?? 
      gbt.est <- c(fit$beta[,1],0)
      cor.r[ii,k] <- cor(gbt.est, lambda.vec, method = 'kendall')
      
    }
    # end of max.k iter
    #cat( which.max( colMeans(cor.r[1:ii, , drop = F], na.rm = T)), '\n')
    #cat(which.max( colMeans(cor.topo[1:ii,, drop = F], na.rm = T)), '\n')
    #cat( mean(cor.naive[1:ii]), '\n')    
  }
  cor.naive_list[[tn_i]] = cor.naive
  cor.r_list[[tn_i]] = cor.r
}

cor.r_list5 = cor.r_list
cor.naive_list5 = cor.naive_list



colMeans(cor.r_list1[[1]] , na.rm=T)
colMeans(cor.r_list1[[2]] , na.rm=T)
colMeans(cor.r_list1[[3]] , na.rm=T)

colMeans(cor.r_list3[[1]] , na.rm=T)
colMeans(cor.r_list3[[2]] , na.rm=T)
colMeans(cor.r_list3[[3]] , na.rm=T)

colMeans(cor.r_list5[[1]] , na.rm=T)
colMeans(cor.r_list5[[2]] , na.rm=T)
colMeans(cor.r_list5[[3]] , na.rm=T)

setwd("E:\\rank consistency\\simulation\\Rdata\\")
save.image('.\\Simulation2_dh_170613.rdata')    

## plot 그리기 
setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim5_dh_ver2_1.jpg')
plot(0:9 , colMeans(cor.r_list1[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
lines(0:9 , colMeans(cor.r_list1[[2]] , na.rm=T) , type='b' , lty=2)
lines(0:9 , colMeans(cor.r_list1[[3]] , na.rm=T) , type='b' , lty=3)
abline(v=1 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver2_3.jpg')
plot(0:9 , colMeans(cor.r_list3[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
lines(0:9 , colMeans(cor.r_list3[[2]] , na.rm=T) , type='b' , lty=2)
lines(0:9 , colMeans(cor.r_list3[[3]] , na.rm=T) , type='b' , lty=3)
abline(v=3 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver2_5.jpg')
plot(0:9 , colMeans(cor.r_list5[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
lines(0:9 , colMeans(cor.r_list5[[2]] , na.rm=T) , type='b' , lty=2)
lines(0:9 , colMeans(cor.r_list5[[3]] , na.rm=T) , type='b' , lty=3)
abline(v=5 , col='red')
dev.off()


################################################################################
## 하나씩 따로 그리기
setwd("E:\\rank consistency\\simulation\\Resdata\\")
jpeg('.\\sim5_dh_ver3_1_500.jpg')
plot(0:9 , colMeans(cor.r_list1[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=1 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_1_5000.jpg')
plot(0:9 , colMeans(cor.r_list1[[2]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=1 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_1_50000.jpg')
plot(0:9 , colMeans(cor.r_list1[[3]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=1 , col='red')
dev.off()


jpeg('.\\sim5_dh_ver3_3_500.jpg')
plot(0:9 , colMeans(cor.r_list3[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation' , 
     cex.lab=1.8, cex.axis=1.5, cex.sub=2.3 , lwd=1.5)
abline(v=3 , col='red' , lwd=1.5)
dev.off()

jpeg('.\\sim5_dh_ver3_3_5000.jpg')
plot(0:9 , colMeans(cor.r_list3[[2]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation' , 
     cex.lab=1.8, cex.axis=1.5, cex.sub=2.3 , lwd=1.5)
abline(v=3 , col='red' , lwd=1.5)
dev.off()

jpeg('.\\sim5_dh_ver3_3_50000.jpg')
plot(0:9 , colMeans(cor.r_list3[[3]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation' , 
     cex.lab=1.8, cex.axis=1.5, cex.sub=2.3 , lwd=1.5)
abline(v=3 , col='red' , lwd=1.5)
dev.off()


jpeg('.\\sim5_dh_ver3_5_500.jpg')
plot(0:9 , colMeans(cor.r_list5[[1]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=5 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_5_5000.jpg')
plot(0:9 , colMeans(cor.r_list5[[2]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=5 , col='red')
dev.off()

jpeg('.\\sim5_dh_ver3_5_50000.jpg')
plot(0:9 , colMeans(cor.r_list5[[3]] , na.rm=T) , ylim=c(0.7,1) , type='b' , lty=1 , 
     xlab='threshold level' , ylab='correlation')
abline(v=5 , col='red')
dev.off()

