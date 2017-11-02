# parallel_build_mim

This functions build a pairwise correlation matrix or a mutual information matrix given a dataset.
If the dataset is an NxM matrix it builds an NxN matrix
The function uses the foreach and doMC function for parallelization purpuse
Mutual information is evaluated by using the infotheo package
