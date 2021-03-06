// [[Rcpp::interfaces(r, cpp)]]
#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


bool isNonIncreasing(arma::vec x){
  /*Function to check a sequence is if nonincreasing or not
   * input sequence x as a vector
   * returns true: if nonincreasing
   *         false: otherwises
   */
  bool flag = false;
  for(int i=0;i+1<x.n_elem;i++){
    if(x(i)>=x(i+1))
      flag = true;
    else
      flag = false;
  }
  return flag;
}
//[[Rcpp::export]]
arma::vec FastProxSL1(arma::vec y, arma::vec lambda){
  /*Function calculates prox operator using ALgorithm3
   * input: vector y - Nonnegative and nonincreasing 
   *        vector lambda - Nonnegative and nonincreasing 
   * returns: vector x
   */
  
  arma::vec x = y - lambda;
  int n = x.n_elem;
  int i, j, k;
  double avg_y, avg_lambda;
  
  while(!isNonIncreasing(x)){
    for(i=n-2;i>=0;i--){
      j = i+1;
      
      while(j<n && x(j-1)<x(j)) //Identify nondecreasing and nonconstant subsequences
        j++;
      
      j--;
      
      if(i == j)
        continue;
      
      //Calculate with averages
      avg_y =0;
      avg_lambda =0;
      
      for(k =i;k<j;k++){
        avg_y += y(k);
        avg_lambda += lambda(k);
      } 
      
      avg_y = avg_y /(j-i+1);
      avg_lambda = avg_lambda /(j-i+1);
      
      //Replace with averages  
      for(k =i;k<j;k++){
        y(k) = avg_y;
        lambda(k) = avg_lambda;
      }
      //Update x
      x = y - lambda;
    }
  }
  for(i=0;i<n;i++)
    if(x(i)<0)
      x(i)=0;
    return x;
}