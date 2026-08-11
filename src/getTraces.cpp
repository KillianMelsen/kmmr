#include <Rcpp.h>
#include <algorithm> // Required for std::lower_bound
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



// Helper function to safely fetch a value from a sparse matrix
// Now completely agnostic to dsCMatrix ("L" or "U") and dgCMatrix!
inline double get_sparse_val(int r, int c, const int* p, const int* i, const double* x) {
  
  // 1. First, try looking for the exact (r, c) coordinate
  int start = p[c];
  int end = p[c + 1];
  
  if (start < end) {
    const int* pos = std::lower_bound(i + start, i + end, r);
    if (pos != (i + end) && *pos == r) {
      return x[pos - i];
    }
  }
  
  // 2. If not found, and since the matrix is symmetric, try looking up (c, r)
  // This safely finds the missing triangle in a dsCMatrix without knowing uplo.
  if (r != c) {
    start = p[r];
    end = p[r + 1];
    
    if (start < end) {
      const int* pos = std::lower_bound(i + start, i + end, c);
      if (pos != (i + end) && *pos == c) {
        return x[pos - i];
      }
    }
  }
  
  // 3. If neither is found, it is a structural zero
  return 0.0;
}

// [[Rcpp::export]]
NumericMatrix getTracesSparse(IntegerVector C_p, IntegerVector C_i, NumericVector C_x, 
                              IntegerVector G_p, IntegerVector G_i, NumericVector G_x, 
                              int p, int d) {
  
  NumericMatrix traces(p, p);
  
  // Extract raw C++ pointers from Rcpp vectors for maximum speed
  const int* G_p_ptr = INTEGER(G_p);
  const int* G_i_ptr = INTEGER(G_i);
  const double* G_x_ptr = REAL(G_x);
  
  for (int col = 0; col < C_p.length() - 1; ++col) {
    int block_j = col / d;
    int in_j = col % d;
    
    int start = C_p[col];
    int end = C_p[col + 1];
    
    for (int idx = start; idx < end; ++idx) {
      int row = C_i[idx];
      int block_i = row / d;
      int in_i = row % d;
      
      // Fast binary-search lookup in sparse G (now uplo-agnostic)
      double g_val = get_sparse_val(in_i, in_j, G_p_ptr, G_i_ptr, G_x_ptr);
      
      // Skip the multiplication if G is structurally zero here
      if (g_val == 0.0) continue;
      
      double val = C_x[idx] * g_val;
      
      // C_sub is still assumed to be lower triangular (from your C matrix)
      if (block_i == block_j) {
        if (in_i > in_j) {
          val *= 2.0; 
        }
      }
      
      traces(block_j, block_i) += val;
    }
  }
  return traces;
}


