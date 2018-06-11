## description
## investigate the variance of the proposed estimator
rm(list = ls()); gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
require('glmnet')
library(igraph)
library(ggplot2)
library(dplyr)
source('./lib/sim.R')
################################################################################
max.k = 10
p = 10
kn <- 7  ## d
rn <- 3   ## n_s
df = 1    ## degree of freedom
counter = 1
sim.iter = 200
source('./lib/exe-2.R')
tn_vec = c(500,5000,50000)
cor.naive_list = list()
cor.r_list = list()
k_fold = 5
tn_i = 3

# simulation: number of obs.
  tn = tn_vec[tn_i]  ## tn 정의 (전체 rank pair의 수.)
  cat ('total n:' , tn , '\n')
  cor.naive<- rep(0,sim.iter) ## BT를 이용한 kendall's tau 저장하는 벡터 
  
  
  gbt_trueMat <- matrix(NA,sim.iter,max.k)
  gbt_truevec <- rep(NA, sim.iter)
  bt_truevec <- rep(NA, sim.iter)
  
  cor.cv <- matrix(0,sim.iter,max.k)
  cor.cv.list = vector(mode = 'list', length = sim.iter)
  ii = 1
  for  (ii in 1:sim.iter)
  {
    cat(' ',ii,'-th iteration\n')
    set.seed(ii+123) ## set seed
    Qmat = sparse_gen_fun(dmat, kn, rn, tn)
    gen_fit = gen_sim_fun(Gmat, Qmat)
    cvec<-(0:(max.k-1))/tn  ## cvec : threshold cval
    
    # cross validation
    fit = cv.gbt_fun(gen_fit, cvec, k_fold, lambda.vec)
    cv_mean_vec = colMeans(fit, na.rm = TRUE)
    cor.cv.list[[ii]] = fit
    cor.cv[ii,] =  cv_mean_vec
    
    # cor.r
    for (k in 1:length(cvec))
    {
      cval = cvec[k]
      gbt_fit = gbt_fun(gen_fit, cval, lambda.vec)
      gbt_trueMat[ii,k] = gbt_fit$cor
    }
    
    # cv selection
    min_idx <-which.max(cv_mean_vec)
    cval <-  cvec[min_idx]

    gbt_fit = gbt_fun(gen_fit, cval, lambda.vec)
    bt_fit = bt_fun(gen_fit, lambda.vec)
    
    gbt_truevec[ii] = gbt_fit$cor
    bt_truevec[ii] = bt_fit$cor
  } 
boxplot(bt_truevec, gbt_truevec, names = c("BT", "gBT2"),
        col = 'lightblue')      
#load("sim_result_5-3.rdata")      
# fig 1
labsize = 2.2
load("simulation-5-sc.rdata")
cor.cv_jj = cor.cv
load("sim_result_5-1.rdata")
png(filename = "sim5_1_500_JJ.png", width = 600, height = 600)
boxplot(cor.cv_jj[,1], bt_truevec, gbt_trueMat[,1], names = c("SC","BT", "gBT2"),
        col = 'lightblue', ylim = c(0.60,1),
        ylab = "correlation", cex.lab = labsize,
        cex.axis = labsize)      
dev.off()

png(filename = "sim5_1_500_CV.png", width = 600, height = 600)
boxplot(gbt_trueMat[,1], gbt_truevec, names = c("gBT2", "gBT2-CV"),
        col = 'lightblue', ylim = c(0.60,1),
        ylab = "correlation", cex.lab = labsize,
        cex.axis = labsize)      
dev.off()

# fig 2
load("simulation-5-sc.rdata")
load("sim_result_5-2.rdata")      
png(filename = "sim5_1_5000_JJ.png", width = 600, height = 600)
boxplot(cor.cv_jj[,2], bt_truevec, gbt_trueMat[,1], names = c("SC","BT", "gBT2"),
        col = 'lightblue', ylim = c(0.60,1),
        ylab = "correlation", cex.lab = labsize,
        cex.axis = labsize)      

dev.off()

png(filename = "sim5_1_5000_CV.png", width = 600, height = 600)
boxplot(gbt_trueMat[,1], gbt_truevec, names = c("gBT2", "gBT2-CV"),
        col = 'lightblue', ylim = c(0.60,1),
        ylab = "correlation", cex.lab = labsize,
        cex.axis = labsize)      
dev.off()

# fig 3
load("simulation-5-sc.rdata")
load("sim_result_5-3.rdata")      
png(filename = "sim5_1_50000_JJ.png", width = 600, height = 600)
boxplot(cor.cv_jj[,3], bt_truevec, gbt_trueMat[,3], names = c("SC","BT", "gBT2"),
        col = 'lightblue', ylim = c(0.60,1),
        ylab = "correlation", cex.lab = labsize,
        cex.axis = labsize)      
dev.off()

png(filename = "sim5_1_50000_CV.png", width = 600, height = 600)
boxplot(gbt_trueMat[,1], gbt_truevec, names = c("gBT2", "gBT2-CV"),
        col = 'lightblue', ylim = c(0.60,1),
        ylab = "correlation", cex.lab = labsize,
        cex.axis = labsize)      
dev.off()


