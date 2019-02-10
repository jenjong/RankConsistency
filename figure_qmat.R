

plot(x = 0,y = 0, 
     xlim = c(0,1), ylim = c(0,1), 
     type ='n', xlab = 'x', ylab = 'y')
abline(v = 0.1, lty = 2)
abline(h = 0.1, lty = 2)
points(0.1, 0.1)
p = 10
ddirichlet<- function(x,alpha.vec)
{
  xd = 1-sum(x)
  v<- sum((alpha.vec[1:2]-1)*log(x)) + (alpha.vec[3]-1)*log(xd)
  v1<- gamma(sum(alpha.vec))/prod(gamma(alpha.vec))*exp(v)
  v1  
}

  alpha.vec1 = c(10, 8, 3)
  alpha.vec2 = c(3, 3, 10)
  alpha.vec3 = c(10, 3, 4)
  pi.var1 = 1/3
  pi.var2 = 1/3
  pi.var3 = 1/3
  
  z = seq(0, 1, length = p+2)
  z = z[2: ((length(z)-1))]
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
  #install.packages(rgl)
  library(rgl)
  persp3d(z,z,dmat, col='lightblue', xlab = expression(x[1]), 
          ylab = expression(x[2]), zlab = expression(f(x[1],x[2])),
          cex.axis = 1.5)
  
  dmat <- dmat/sum(dmat)
  image(dmat, col = gray.colors(10, start = 0, end = 1))
  dmat2 <- matrix(0,p,p)
  for (j in 1:p) dmat2[,j] <- rev(dmat[j,])
  dmat2 <- dmat2 + t(dmat2) 
  Qmat = dmat2
  
  imat <- matrix(0,p,p)
  for (j in 1:p)
    for (k in 1:p)
      imat[k,p+1-j]<- Qmat[j,k]
  a<- max(imat)
  for (i in 1:9)
  {
    for (j in 1:(10-i))
    {
      imat[i,j]<- a
    }
  }
  image(x = 1:p, y = 1:p, imat, col = gray.colors(100,start = 0, end = 1),
        main = 'k', ylab = 'j', xlab = "", axes=FALSE, font.lab = 2, 
        cex.lab = 1.1)
  axis(2, at=p:1,labels=1:p,
       col.axis="blue", las=2, cex.axis=0.9, tck=-.01)  
  axis(3, at=1:p,labels=1:p,
       col.axis="blue", las=2, cex.axis=0.9, tck=-.01)
  