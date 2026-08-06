#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getTraces(IntegerVector p_ptr, IntegerVector i_ptr, NumericVector x_ptr, NumericMatrix G, int p, int d) {
  
  NumericMatrix traces(p, p);
  
  for (int col = 0; col < p_ptr.length() - 1; ++col) {
    int block_j = col / d;
    int in_j = col % d;
    
    int start = p_ptr[col];
    int end = p_ptr[col + 1];
    
    for (int idx = start; idx < end; ++idx) {
      int row = i_ptr[idx];
      int block_i = row / d;
      int in_i = row % d;
      
      double val = x_ptr[idx] * G(in_i, in_j);
      
      // If we are calculating a block on the diagonal...
      if (block_i == block_j) {
        // ... double the lower off-diagonal elements to account for the missing upper half
        if (in_i > in_j) {
          val *= 2.0; 
        }
      }
      
      // Assign to (block_j, block_i) to form the upper triangle
      traces(block_j, block_i) += val;
    }
  }
  return traces;
}


