#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin,int a,int b,int n) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1,n,y)[0];
      y = rbeta(1, x + a, n-x+b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
