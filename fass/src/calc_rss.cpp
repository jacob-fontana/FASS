#include <Rcpp.h>
using namespace Rcpp;

//' Fast Residual Sum of Squares Calculation
//'
//' Efficiently computes residual sum of squares across active covariates.
//' @param X Design matrix
//' @param y Target response vector (or centered/intercept-adjusted response)
//' @param beta_act Active regression coefficients
//' @param idx_act 1-based indices of active features in X
//' @return Sum of squared residuals
// [[Rcpp::export]]
double calc_rss_cpp(NumericMatrix X, NumericVector y, NumericVector beta_act, IntegerVector idx_act) {
  int n = X.nrow();
  int p = X.ncol();
  int k = idx_act.size();
  
  // Allocate a single vector for the residuals
  NumericVector resid = clone(y); 
  
  // Cache-friendly column-major loop
  for(int j = 0; j < k; ++j) {
    int col = idx_act[j] - 1; // C++ is 0-indexed, R is 1-indexed
    if (col < 0 || col >= p) continue; // Bounds check protection
    double b = beta_act[j];
    
    for(int i = 0; i < n; ++i) {
      resid[i] -= X(i, col) * b;
    }
  }
  
  // Calculate the sum of squares with floating-point stability
  double rss = 0.0;
  for(int i = 0; i < n; ++i) {
    rss += resid[i] * resid[i];
  }
  
  return rss;
}
