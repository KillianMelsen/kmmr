#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
void inplaceUpdateSIS(S4 Ciss, S4 tW_mat, S4 tQ_mat) {
  // 1. Extract sparse slots directly from the dsCMatrix
  NumericVector x_ptr = Ciss.slot("x");
  IntegerVector p_ptr = Ciss.slot("p");
  IntegerVector i_ptr = Ciss.slot("i");
  
  // 2. Map S4 dense matrices to Armadillo (zero-copy mapping)
  // tW_mat is (k x N_s)
  NumericVector tW_x = tW_mat.slot("x");
  IntegerVector tW_dim = tW_mat.slot("Dim");
  arma::mat tW(tW_x.begin(), tW_dim[0], tW_dim[1], false, true);
  
  // tQ_mat is (k x N_s)
  NumericVector tQ_x = tQ_mat.slot("x");
  IntegerVector tQ_dim = tQ_mat.slot("Dim");
  arma::mat tQ(tQ_x.begin(), tQ_dim[0], tQ_dim[1], false, true);
  
  int ncol = p_ptr.length() - 1;
  
  // 3. Iterate and update in-place using contiguous memory
  for (int j = 0; j < ncol; ++j) {
    int start = p_ptr[j];
    int end = p_ptr[j + 1];
    
    for (int idx = start; idx < end; ++idx) {
      int row = i_ptr[idx]; // Corresponds to index i
      
      // Blazing fast vectorized dot product down the columns
      x_ptr[idx] += arma::dot(tW.col(row), tQ.col(j));
    }
  }
}


