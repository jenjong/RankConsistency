
########################################
# glmnet design matrix and reponse matrix
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

###############################  
# Gmat : population version 
lambda.vec = seq(10,0,length = p)
#lambda.vec= qt(seq(0.9, 0.5, length = p), df = 1)
Gmat = pt( outer(lambda.vec, lambda.vec, FUN = '-'), df = df)
##############################

##########################################################################################

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
dmat <- dmat/sum(dmat)
