#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

//' Vectorized Continuous Ranked Probability Score and Quantile Metrics
//'
//' Computes sample mean, quantiles (2.5%, 97.5%), and empirical CRPS.
//' @param preds An n x M matrix of posterior predictive draws
//' @param y Target response vector of length n
//' @return List containing crps, mean, q2_5, and q97_5 vectors
// [[Rcpp::export]]
List compute_gibbs_metrics_cpp(NumericMatrix preds, NumericVector y) {
  int n = preds.nrow();
  int M = preds.ncol();
  
  NumericVector crps(n);
  NumericVector p_mean(n);
  NumericVector p_q25(n);
  NumericVector p_q975(n);
  
  int idx_25 = std::max(0, (int)(0.025 * M));
  int idx_975 = std::min(M - 1, (int)(0.975 * M));
  
  for (int i = 0; i < n; i++) {
    std::vector<double> x(M);
    double sum_x = 0;
    for (int j = 0; j < M; j++) {
      x[j] = preds(i, j);
      sum_x += x[j];
    }
    p_mean[i] = sum_x / M;
    
    std::sort(x.begin(), x.end());
    
    p_q25[i] = x[idx_25];
    p_q975[i] = x[idx_975];
    
    double term1 = 0;
    double term2 = 0;
    for (int j = 0; j < M; j++) {
      term1 += std::abs(x[j] - y[i]);
      term2 += (2.0 * (j + 1) - 1.0 - M) * x[j];
    }
    term1 /= M;
    term2 /= ((double)M * M);
    crps[i] = term1 - term2;
  }
  
  return List::create(Named("crps") = crps,
                      Named("mean") = p_mean,
                      Named("q2_5") = p_q25,
                      Named("q97_5") = p_q975);
}
