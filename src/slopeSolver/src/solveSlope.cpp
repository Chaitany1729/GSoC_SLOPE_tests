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
//[[Rcpp::export]]

arma::vec solveSlope(arma::mat x, arma::vec y, arma::vec lambda){

  int p = x.n_cols;
  int iterations = 10000;
  double rho = 10;
  
  arma::vec beta(p), z, w;
  arma::mat I(p,p);
  I.eye();
  beta.fill(0);
  z = beta;
  w = beta;
  

  arma::mat y_hat = x.t()*y;
  arma::mat t = inv(x.t()*x + rho*I);
  
  for(int i=0;i<iterations;i++){
    
    beta = t*(y_hat+rho*(z-w));
    z = beta + w;
    arma::vec z_sign = sign(z);
    z = abs(z);
    arma::uvec z_index = sort_index(z,"descend");
    z = (z(z_index)).eval();
    
    z = FastProxSL1(z, (1/rho)*lambda);
    z(z_index) = z;
    z %= z_sign;
    
    w = w + beta - z;
  }
  
  return beta;
}
