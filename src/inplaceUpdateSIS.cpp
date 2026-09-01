#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void inplaceUpdateSIS(NumericVector x_ptr, IntegerVector p_ptr, IntegerVector i_ptr,
                      NumericMatrix W_mat, NumericMatrix Q) {
  int ncol = p_ptr.length() - 1; // Number of columns pf Css, equal to the number of sparse equations
  int k = W_mat.ncol(); // Number of columns of Q, equal to the number of dense equations
  
  // Iterate over columns:
  for (int j = 0; j < ncol; ++j) {
    int start = p_ptr[j];
    int end = p_ptr[j + 1];
    
    // Iterate through structural non-zeros in column j:
    for (int idx = start; idx < end; ++idx) {
      int i = i_ptr[idx]; // Row index
      
      // Compute dot product of row i of W_mat and row j of Q.
      // This is equal to computing element [i, j] of Q %*% solve(SC) %*% t(Q):
      double update_val = 0.0;
      for (int m = 0; m < k; ++m) {
        update_val += W_mat(i, m) * Q(j, m);
      }
      
      // Add directly to the existing SIS value in-place
      x_ptr[idx] += update_val;
    }
  }
}