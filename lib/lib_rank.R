IWLS.ridge <- function(X,true.prob,pij,accuracy=10e-10,maxitration=2000)
{
  p = ncol(X)
  tX = t(X)
  initial=rep(0,p)
  pv=1/(1+exp(-X%*%initial))
  derive.1=tX%*%((true.prob-pv)*pij)
  derive.2=-tX%*%diag(drop(pv)*(1-drop(pv))*pij)%*%X
  beta.old=initial-solve(derive.2)%*%derive.1
  bool <- T
  beta.new <- rep(0,p)
  Step <- 1
  Beta <- drop(beta.old)
  maxitr=1
  while(bool){
    #cat(Step, "th : ", beta.old, "\n")
    pv=1/(1+exp(-X%*%beta.old))
    derive.1=tX%*%((true.prob-pv)*pij)
    derive.2=- tX%*%diag(drop(pv)*(1-drop(pv))*pij)%*%X
    beta.new=beta.old-solve(derive.2)%*%derive.1
    Beta=c(Beta,drop(beta.new))
    maxitr=maxitr+1
    if (max(abs(drop(beta.new)-drop(beta.old)))<accuracy||maxitr==maxitration) bool=F
    beta.old=beta.new
  }
  beta.old
}

desg_fun = function(pl_num)
{
  X = matrix(0, choose(pl_num,2), pl_num )

  k = 1
  for ( i in ( 1 :(pl_num-1) ) )
  {
    for ( j in (i+1):(pl_num) )
    {
      X[k,i] = 1
      X[k,j] = -1
      k = k + 1
    }

  }
  X = X[,-pl_num]
  return(X)
}

mix.fun = function(n, pi = 0.5, mu1 = 0, mu2 = 1, sd1 = 1, sd2 = 0.2)
{
  x1 = rnorm(n,mu1,sd1)
  x2 =  rnorm(n,mu2,sd2)
  x3 = runif(n)
  x1[x3>pi] = x2[x3>pi]
  x1
}

### 
gen.mat= function(x, type = 'probability')
    {
        p =  (1 + sqrt(1+8*length(x)) )/2
        a = matrix(0,p,p)
        idx = 1
        for( i in 1:(p-1))
        {
            for (j in (i+1):p)
            {
            a[i,j]= x[idx]
            idx = idx + 1
            }
        }
        if (type == 'probability')
        {
               for (i in 2:p)
               {
                    for ( j in 1:(i-1) )
                    {
                    a[i,j] = 1 - a[j,i]
                    }
               }
        }

        if (type == 'number')
        {
                a = a + t(a)
        }

    return(a)
    }

rank.fun = function(x, p)
    {
    y = p - t( apply(x,1,rank)) + 1
    return(y)
    }

bal.mat = function(x, j1,j2)
        {
        idx = (x[j1,] != 0 &  x[j2,] != 0)
        if ( sum (idx)>0 ) 
        {
        x[j1,!idx] = x[j2,!idx] = x[!idx,j1] = x[!idx,j2]=0
        tmp = (x[j1,idx] + x[j2,idx])/2
        x[j1,idx] <- x[j2,idx] <- x[idx,j1] <- x[idx,j2] <- tmp
        for ( i in 1:ncol(x)) x[i,i] = 0
        }
        return(x)
        }

gen.vec = function(x,p)
        {
        y = c()
        idx = 1
        for (i in 1:(p-1))
        {
            for (j in (i+1):p)
            {
            y[idx] = x[i,j]
            idx = idx + 1
            }
        }
        return(y)
        }
### SVM for perfect separation ###
opt.sep = function(x,y,p)
        {
        y = 2*y - 1
        Dmat = diag(1,p-1)
        dvec = rep(0,p-1)
        Amat = t(X*y)
        bvec = rep(1,nrow(x))
        fit = solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
        beta.vec = c(fit$solution,0)
        return(beta.vec)
        }

####
pIWLS = function(X,true.prob,pij,accuracy=10e-10,maxitration=200)
        {
        p = ncol(X) + 1
        pij.mat = gen.mat(pij,type='number')
        my = c()
        idx = 1
        for (i1 in 1:(p-1))
        {
            for (i2 in (i1+1):p )
            {
            mij.mat = bal.mat(pij.mat,i1,i2)
            mij = gen.vec(mij.mat,p)
            beta.vec= IWLS.ridge(X,true.prob,mij,accuracy=accuracy,maxitration=maxitration)
            beta.vec = c(beta.vec,0)
            my[idx] = as.integer(beta.vec[i1] - beta.vec[i2]>0)
            idx = idx + 1
            }
        }
        beta.vec = opt.sep(X,my,p)
        return(beta.vec)
        }
