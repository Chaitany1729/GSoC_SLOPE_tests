library('Rcpp')
library('SLOPE')

sourceCpp('FastProxSL1.cpp')

length_list = as.integer(runif(1, 100, 10000))
time_fastprox = c()
time_proxsort = c()
for (vector_length in length_list){
  y = runif(vector_length, 0, 1000)
  y = sort(y, decreasing = TRUE)
  lambda = runif(vector_length, 0, 1000)
  lambda = sort(lambda, decreasing = TRUE)
  t = proc.time()['elapsed']
  FastProxSL1(y, lambda)
  t = proc.time()['elapsed'] - t
  time_fastprox = c(time_fastprox, t)
  t = proc.time()['elapsed']
  prox_sorted_L1(y, lambda)
  t = proc.time()['elapsed'] - t
  time_proxsort = c(time_proxsort, t)
}



