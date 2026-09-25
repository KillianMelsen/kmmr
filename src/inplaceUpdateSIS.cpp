#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void inplaceUpdateSIS(NumericVector x_ptr, IntegerVector p_ptr, IntegerVector i_ptr, NumericMatrix W_mat, NumericMatrix Q) {
  int ncol = p_ptr.length() - 1; // Number of columns of Ciss, equal to the number of sparse equations
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

// [[Rcpp::export]]
void inplaceUpdateSISnew(NumericVector x_ptr, IntegerVector p_ptr, IntegerVector i_ptr, NumericMatrix tW_mat, NumericMatrix tQ) {
  int ncol = p_ptr.length() - 1; 
  // k is now the number of rows of the transposed matrices
  int k = tW_mat.nrow(); 
  
  // Iterate over columns of Css:
  for (int j = 0; j < ncol; ++j) {
    int start = p_ptr[j];
    int end = p_ptr[j + 1];
    
    // Grab a contiguous memory pointer to the j-th column of tQ
    const double* q_col = &tQ(0, j);
    
    // Iterate through structural non-zeros in column j:
    for (int idx = start; idx < end; ++idx) {
      int i = i_ptr[idx]; 
      
      // Grab a contiguous memory pointer to the i-th column of tW_mat
      const double* w_col = &tW_mat(0, i);
      
      double update_val = 0.0;
      
      // Blazing fast contiguous memory access
      for (int m = 0; m < k; ++m) {
        update_val += w_col[m] * q_col[m];
      }
      
      x_ptr[idx] += update_val;
    }
  }
}