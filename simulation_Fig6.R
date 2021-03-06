## description
rm(list = ls()); gc()
setwd("C:/Users/Jeon/Documents/GitHub/RankConsistency")
require('glmnet')
library(igraph)
library(ggplot2)
source('./lib/car_lib.R')
source('./lib/lib_rank.R')
source('./lib/sim.R')
source('./lib/real_lib.R')
max.k = 10
p = 10
kn <- 7  ## d
rn <- 3   ## n_s
df = 1    ## degree of freedom
counter = 1
source('./lib/exe-2.R')
# here we already select the two pairs
dmat[c(11, 81)] = 0
tn = 50000
cor.naive_list = list()
cor.r_list = list()

sim.iter = 50
cor.bt = cor.gbt = cor.sr = c()
corMat.gbt = NULL
cut_vec = 1:8
jj = 1
for (jj in 1:length(cut_vec))
{
  cat(jj,'\n')
  cut_var = cut_vec[jj]
  ii = 1
  cor.gbt = c()
  for ( ii in 1:sim.iter)
  {
    set.seed(ii+1105)
    #cat(ii,'\n')
    ### generating number of comparison using dmat
    ### ouput: Qmat
    dmat1 <- dmat
    u.idx <- which( dmat > 0)
    sel.u.idx<- sample(u.idx, kn-2)
    dmat1[sel.u.idx]  <- 0
    dmat1 <- dmat1/sum(dmat1)
    d.sample <- drop (rmultinom(1, tn-rn*kn, prob = c(dmat1)))
    d.sample[sel.u.idx] <- rn
    dmat1 <- matrix(d.sample, p , p)
    Qmat <- matrix(0, p, p )
    for (j in 1:p) Qmat[,j] <- rev(dmat1[j,])
    Qmat <- Qmat + t(Qmat)
    ##############################
    # generating the result of win-loss ratio using Gmat(true probability matrix) on Qmat (number of matching)
    # output: Gmat.hat
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
    w_mat  = Gmat.hat
    Gmat.hat <- Gmat.hat/Qmat
    Gmat.hat[!is.finite(Gmat.hat)] = 0
    ###############################
    # Generating Qpmat that consists of q_{jk}
    n = sum(Qmat)
    Qpmat = Qmat/n*2
    
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
    
    Qmat_fit <- list()
    Qmat_fit$x = x
    Qmat_fit$y = y
    Qmat_fit$Qpmat = Qpmat
    Qmat_fit$Gmat_hat = Gmat.hat
    Qmat_fit$n = n
    Qmat_fit$sc_alarm = FALSE
    Qmat_fit$Qmat = Qmat
    Qmat_fit$Wmat = w_mat
    
    
    #bt_est <- btFun(Qmat_fit)
    #sr_est <- srFun(Qmat_fit)
    gbt_fit <- gbtFun(Qmat_fit, cut_v = cut_var/tn, 'balance')
    gbt_est = gbt_fit$gbt_est
    if (is.null(gbt_est)) 
    {
      cor.gbt[ii] = NA
      next
    }
    cor.gbt[ii] = cor(gbt_est, lambda.vec, method = 'kendall')
  }
  corMat.gbt = cbind(corMat.gbt, cor.gbt)
}
save.image("fig6_0805_3.rdata")  

exp.par = 3
png('./fig6_0805_3.PNG', width = 960, height = 960)  
par( mai = c(4,4,4,2)/3)  
plot(0:7, colMeans(corMat.gbt, na.rm = T), 
     ylim = c(0.7,1),
     type = 'b',
     xlab  = 'threshold level',
     ylab = 'correlation', 
     cex.lab = exp.par,
     cex.axis = exp.par,
     lwd = 2,
     cex = 2
    )
abline(v = 3, col = "red", lwd = 2)
dev.off()
# 
# png('./fig6_0805_3.PNG', width = 960, height = 960)  
# par( mai = c(4,4,4,2)/3)  
# exp.par = 3
#   boxplot(cor.bt, cor.sr, cor.gbt, ylim = c(0.6,1),
#           names = c("BT", "SC", "gBT2"), 
#           col = 'lightblue', 
#           ylab = 'correlation', 
#           cex.lab = exp.par,
#           cex.axis = exp.par)
# dev.off()
