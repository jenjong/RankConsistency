
genFun = function(mat.xinfo, mat.yinfo,p)
{ 
  x = NULL
  y = c()
  nij = matrix(0,p,p)
  list.vec = NULL
  idx = 1
  for (i in 1:length(mat.yinfo))
  {
    cat(i,'\n')
    tmpx<-mat.xinfo[i,1:mat.yinfo[i]]
    
    for (j in 1:(length(tmpx)-1))
    {
      for (k in j:length(tmpx) )
      {
        if (tmpx[j] == tmpx[k]) next
        
        if (tmpx[j]<tmpx[k])
        {
          y<-c(y,1)
          xvec = rep(0,p)
          xvec[tmpx[j]] = 1 ; xvec[tmpx[k]] = -1 
          x <- rbind(x,xvec)
          nij[tmpx[j],tmpx[k]] = nij[tmpx[j],tmpx[k]] + 1
        }
        
        if (tmpx[j]>tmpx[k])
        {
          y<-c(y,0)
          xvec = rep(0,p)
          xvec[tmpx[j]] = -1 ; xvec[tmpx[k]] = 1
          x <- rbind(x,xvec)
          nij[tmpx[k],tmpx[j]] = nij[tmpx[k],tmpx[j]] + 1
        }
      }
    }
    if (i%%1000 == 0 )  {
      list.vec[[idx]] = x
      x = NULL
      idx = idx + 1
    }
  }
  list.vec[[idx+1]]  = x
  Nij = nij + t(nij)
  xmat = NULL
  for ( i in 1:length(list.vec)) xmat = rbind(xmat, list.vec[[i]])
  return(list(xmat = xmat, yvec = y, nij = Nij))
}