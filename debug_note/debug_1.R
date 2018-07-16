# debug: tie -> matrix NA -> 


i = 29
result = sr1.result
method='count'
allowties = F

### code in the function
sc_list <- result
Qmat = Qmat_fit$Qmat
Qpmat = Qmat_fit$Qpmat
Gmat_hat = Qmat_fit$Gmat_hat
p = ncol(Qpmat)
#save.image('debug_1.dbg')
if (method == "count")
{
  mat = matrix(0,p,p);
  for (i in 1:nrow(result))
  {
    j = result[i,1];
    k = result[i,2];
    w_jk = result[i,3];
    if (j == 0) 
    {
      mat[j,k] = 0.5;
      mat[k,j] = 0.5;
      next;
    }
    if (w_jk==1) mat[j,k] = 1  else mat[k,j] = 1;
  }
  
  gbt_est = apply(mat,1,sum);
  gbt_est = gbt_est.copy = rank(gbt_est, ties = 'max')
  gbt_est.unique = unique(gbt_est)
  
  if (!allowties & length(gbt_est.unique) != p)
  {
    gbt_est.table = table(gbt_est)
    gbt_est.table.names= names(gbt_est.table)
    gbt_est.dupRank = 
      as.integer(gbt_est.table.names[which(gbt_est.table>=2)])
    for (i in 1:length(gbt_est.dupRank))
    {
      gbt_est.dupCar = which(gbt_est == gbt_est.dupRank[i])
      # gbt_est.dupCar is always sorted by which function
      gbt_est.dupCar.mat = pair_fun(gbt_est.dupCar)
      for (j in 1:nrow(gbt_est.dupCar.mat))
      {
        j1 = gbt_est.dupCar.mat[j,1]
        j2 = gbt_est.dupCar.mat[j,2]
        idx = result[,1]==j1 & result[,2]==j2
        if (any(idx)) 
        {
          if (result[idx,3] == 0) gbt_est.dupCar.mat[j,] = c(j2, j1)
        } else {
          gbt_est.dupCar.mat[j,] = NA
          next
        }
      }
      
      gbt_est.dupCar.mat = gbt_est.dupCar.mat[!is.na(gbt_est.dupCar.mat[,1]),,drop = FALSE]
      ##
      if (nrow(gbt_est.dupCar.mat) == 0)
      {
        fit = apply(Qmat,1,sum)[gbt_est.dupCar]
        gbt_est.dupCar.irank = gbt_est.dupCar[order(fit, decreasing = T)]
        gbt_est.dupCar.assingRank = 
          (gbt_est.dupRank[i]):(gbt_est.dupRank[i]-length(gbt_est.dupCar)+1)
        gbt_est.copy[gbt_est.dupCar] = 
          gbt_est.dupCar.assingRank[match(gbt_est.dupCar.irank, gbt_est.dupCar)]
        next
      } 
    ##    
      gbt_est.dupCar.mat = matrix(as.character(gbt_est.dupCar.mat),,ncol(gbt_est.dupCar.mat))
      gbt_est.dupCar.dgraph = make_graph(edges = t(gbt_est.dupCar.mat))
      
      # more coding
      if (is_dag(gbt_est.dupCar.dgraph))
      {
        fit = topo_sort(gbt_est.dupCar.dgraph)
        gbt_est.dupCar.irank = as.integer(names(fit))
      } else {
        #fit = apply(Qmat,1,sum)[gbt_est.dupCar]
        fit = apply(Qmat[gbt_est.dupCar,gbt_est.dupCar],1,sum)
        gbt_est.dupCar.irank = gbt_est.dupCar[order(fit, decreasing = T)]
      }
      gbt_est.dupCar.assingRank = 
        (gbt_est.dupRank[i]):(gbt_est.dupRank[i]-length(gbt_est.dupCar)+1)
      gbt_est.copy[gbt_est.dupCar] = 
        gbt_est.dupCar.assingRank[match(gbt_est.dupCar.irank, gbt_est.dupCar)]
    }
  }
  
  return( gbt_est = gbt_est.copy )