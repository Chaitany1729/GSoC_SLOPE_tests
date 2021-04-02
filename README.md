## This repository contains test solutions for R-GSoC ABSOLE project

All the codes are available in src, Results contains all plots and output files.

### Test 1: 
The plots shows the comparative regularization path for three parameters alpha, step and fraction deviance. For higher value of alpha the regularization path dies zero slowly in case of SLOPE than lasso, implying a better control on both FDR and power.

### Test 2:
The package is bulit and checked with `devtools:check()`. The error encountered in the first place in stored in the error file. Those errors were mainly due to packages `caret` and `glmnet` were not installed. Packages installed and the build suceesfully.

### Test 3:
A PR is submitted for the issue 10.

### Test 4:
The algorithm 3: FastProxSL1 in the FastProxSL1.cpp using `RcppArmadillo`. The performance is compared with `SLOPE::prox_sorted_L1()` in the file test4.R (Note: SLOPE::prox_sorted_L1() is deprecated after version 0.1.3 thus older version is required.) 
Results: `prox_sorted_L1()` outperforms `FastProxSL1()`. Considerable difference in execution time is observed for values highe than 10,000. 

### Test 5: 
The apckage in implemented in the directory `slopeSolver`. ADMM algorith is implemented to for soliving purpose. The results are verified witl `SLOPE` package. 
