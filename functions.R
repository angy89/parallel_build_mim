# This functions build a pairwise correlation matrix or a mutual information matrix given a dataset.
# If the dataset is an NxM matrix it builds an NxN matrix
# The function uses the foreach and doMC function for parallelization purpuse
# Mutual information is evaluated by using the infotheo package

par_build_mim = function (dataset, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset)),nCores=nCores) 
{
  require(infotheo)
  require(foreach)
  require(doMC)
  
  if (disc == "equalfreq" || disc == "equalwidth" || disc == 
      "globalequalwidth") 
    dataset <- infotheo::discretize(dataset, disc, nbins)
  if (estimator == "pearson" || estimator == "spearman" || 
      estimator == "kendall") {
    CM = par_cor(exp=dataset,method = estimator, use = "complete.obs",nCores=nCores)
    mim <- CM^2
    diag(mim) <- 0
    maxi <- 0.999999
    mim[which(mim > maxi)] <- maxi
    mim <- -0.5 * log(1 - mim)
  }
  else if (estimator == "mi.mm") 
    estimator = "mm"
  else if (estimator == "mi.empirical") 
    estimator = "emp"
  else if (estimator == "mi.sg") 
    estimator = "sg"
  else if (estimator == "mi.shrink") 
    estimator = "shrink"
  else stop("unknown estimator")
  if (estimator == "mm" || estimator == "emp" || estimator == 
      "sg" || estimator == "shrink") {
    mim <- par_mutinf(exp = dataset, method = estimator,nCores=nCores)
    diag(mim) <- 0
  }
  mim[mim < 0] <- 0
  return(mim)
}



#exp must be discretized
par_mutinf = function(exp,method,nCores = 32){
  nGenes = nrow(exp)
  registerDoMC(nCores)
  
  X = foreach(i=1:nGenes) %dopar%{
    message(i)

    mi = unlist(lapply(X = 1:nGenes,FUN = function(j){
      infotheo::mutinformation(exp[i,],exp[j,],method = method)
    }))
    
    mi
  }
  
  CM = matrix(0,nGenes,nGenes)
  pb = txtProgressBar(min=1,max=nGenes,style=3)
  for(i in 1:nGenes){
    CM[i,] = X[[i]]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  rownames(CM) = colnames(CM) = rownames(exp)
  
  return(CM)
}

par_cor = function(exp,method, use = "complete.obs",nCores = 32){
  nGenes = nrow(exp)
  registerDoMC(nCores)
  
  X = foreach(i=1:nGenes) %dopar%{
    message(i)
    cor(exp[i,],t(exp),method = method,use=use)
  }
  
  CM = matrix(0,nGenes,nGenes)
  pb = txtProgressBar(min=1,max=nGenes,style=3)
  for(i in 1:nGenes){
    CM[i,] = X[[i]]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  rownames(CM) = colnames(CM) = rownames(exp)
  
  return(CM)
}
