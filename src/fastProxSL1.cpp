#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec fastProxSL1(arma::vec y, arma::vec lambda) {
  /*
   This implementation assumes that y and lambda are 
   nonnegative, non increasing sequences provided in vector form.
   
   */
  arma::vec x = y - lambda;
  int n = x.n_elem;
  while(1){
    
    for(int left=n-2;left>=0;left--){
      int right=left+1;
      
      while(right<n && x(right-1)<=x(right)) right++;
      right--;
      
      if(left==right) continue;
      // [left,right] is increasing. Replace with their average
      double new_y = 0;
      for(int i=left;i<=right;i++) new_y += y(i);
      new_y = new_y/(right-left+1);
      for(int i=left;i<=right;i++) y(i) = new_y;
      
      double new_lambda = 0;
      for(int i=left;i<=right;i++) new_lambda += lambda(i);
      new_lambda = new_lambda/(right-left+1);
      for(int i=left;i<=right;i++) lambda(i) = new_lambda;
      
      for(int i=left;i<=right;i++) x(i) = y(i) - lambda(i);
    }
    
    int flag=0;
    
    for(int i=0;i+1<n;i++){
      if( x(i)<x(i+1)) flag=1;
    }
    
    if(flag==0){
      for(int i=0;i<n;i++){
        if(x(i)<0) x(i)=0;
        
      }
      return x;
    }
  }
}