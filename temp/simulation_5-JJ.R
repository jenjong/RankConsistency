## description
## investigate the variance of the proposed estimator
rm(list = ls()); gc()
#setwd("E:\\rank consistency\\simulation\\code\\")
#setwd("C:/Users/uos_stat/Dropbox/A rank-consistency/prog/temp")
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
require('glmnet')
library(igraph)
library(ggplot2)
library(dplyr)
source('./lib/sim.R')

##########################################################################################
max.k = 10
p = 10
kn <- 7  ## d
rn <- 3   ## n_s
df = 1    ## t분포 자유도 
counter = 1
sim.iter = 200    ## 전체 simulation 과정 반복 수 
source('./lib/exe-2.R')
tn_vec = c(500,5000,50000)
cor.naive_list = list()
cor.r_list = list()
k_fold = 5
tn_i = 3
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
    # function: gen_sim_fun, cv_mat_fun, 
    gen_fit <- gen_sim_fun(Gmat, Qmat)
    Gmat.hat_raw = gen_fit$G
    Qmat_raw = geb_fit$Q
    ###
    # k-fold
    cv_mat <-cv_mat_fun(Gmat_obs, Qmat) 
    k_num = 1
    # for (k_num in 1:k_fold)
    # {
    #    
    # }
    
    tmp_te = cv_mat[cv_mat[,k_num] == k_num,]
    tmp_tr = cv_mat[cv_mat[,k_num] != k_num,]
    cv_m = tmp_tr[,1:3]
    cv_table <- cv_table_fun(cv_m)
    Gmat.hat <- cv_table$G
    Qmat <- cv_table$Q
    ## strat cv     
    Gmat.hat <- Gmat.hat/Qmat
    Gmat.hat[!is.finite(Gmat.hat)] = 0
    
    # Generating Qpmat that consists of q_{jk}
    n = sum(Qmat)
    Qpmat = Qmat/n*2
    
    # define the variable to restore the simulation results
    Result = NULL
    Naive.list = list()
    Result.list = list()

    naive.est <- naive_BT_fun(Qpmat, Gmat.hat, p)  
    k = 1
    # Start the algrotihm for each c (threshold variable)
    for (k in 1:length(cvec))
    {
      ######### gBT model ###########
      # set weight-vector
      cval <- cvec[k]
      result <- gbt_step1_fun(Qpmat, Gmat.hat, p, cval)
      Result <- cbind(Result, result[,4])
      Result.list[[k]] <- result
    }
    cor.naive[ii] <- cor(naive.est, lambda.vec, method = 'kendall') 

    
    
    for (k in 1:length(cvec))
    {
      tmp<-Result.list[[k]]
      not0_ind = (tmp[,1]!=0)
      tmp <-tmp[not0_ind, 1:3]
      p.set <-sort(unique(c(tmp[,1:2])))
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

