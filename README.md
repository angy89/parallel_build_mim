This R functions build a pairwise NxN correlation matrix or a mutual information matrix given and NxM dataset

The function uses the foreach and doMC function for parallelization purpuse. Mutual information is computed by using the infotheo package


#Demo
```R
exp = matrix(runif(n = 10000*100,0,1),nrow = 10000,ncol=100)
X = par_build_mim(dataset=exp,nCores = 15, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset))) 
exp = matrix(runif(n = 100*100,0,1),nrow = 100,ncol=100)
X = par_build_mim(dataset=exp,nCores = 15, estimator = "mi.mm", disc = "equalfreq", nbins = sqrt(NROW(exp))) 

```
