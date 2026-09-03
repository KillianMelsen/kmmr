#include <Rcpp.h>
#include <algorithm> // Required for std::lower_bound
using namespace Rcpp;

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

// C is the sparse block of the inverse coefficient matrix.
// K is the dense constant matrix (GRM) passed as a flat vector.
// [[Rcpp::export]]
NumericMatrix getBlockTracesSparseDense(IntegerVector C_p, IntegerVector C_i, NumericVector C_x, NumericVector K_x, int dimVar, int dimNvar) {
  
  NumericMatrix traces(dimVar, dimVar);
  
  // Extract raw pointers for sparse C vectors
  const int* C_p_ptr = INTEGER(C_p);
  const int* C_i_ptr = INTEGER(C_i);
  const double* C_x_ptr = REAL(C_x);
  const double* K_x_ptr = REAL(K_x);
  
  int n_cols = C_p.length() - 1;
  for (int col = 0; col < n_cols; ++col) {
    int block_j = col / dimNvar;
    int in_j = col % dimNvar;
    
    int start = C_p_ptr[col];
    int end = C_p_ptr[col + 1];
    
    for (int idx = start; idx < end; ++idx) {
      int row = C_i_ptr[idx];
      int block_i = row / dimNvar;
      int in_i = row % dimNvar;
      
      // 1D column-major indexing for dense K
      double k_val = K_x_ptr[in_i + in_j * dimNvar];
      double val = C_x_ptr[idx] * k_val;
      
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

// C is the sparse block of the inverse coefficient matrix.
// K is the sparse constant matrix (pedigree based matrix).
// [[Rcpp::export]]
NumericMatrix getBlockTracesSparseSparse(IntegerVector C_p, IntegerVector C_i, NumericVector C_x, IntegerVector K_p, IntegerVector K_i, NumericVector K_x, int dimVar, int dimNvar) {
  
  NumericMatrix traces(dimVar, dimVar);
  
  // Extract raw C++ pointers for both C and K vectors
  const int* C_p_ptr = INTEGER(C_p);
  const int* C_i_ptr = INTEGER(C_i);
  const double* C_x_ptr = REAL(C_x);
  const int* K_p_ptr = INTEGER(K_p);
  const int* K_i_ptr = INTEGER(K_i);
  const double* K_x_ptr = REAL(K_x);
  
  int n_cols = C_p.length() - 1;
  for (int col = 0; col < n_cols; ++col) {
    int block_j = col / dimNvar;
    int in_j = col % dimNvar;
    
    int start = C_p_ptr[col];
    int end = C_p_ptr[col + 1];
    
    for (int idx = start; idx < end; ++idx) {
      int row = C_i_ptr[idx];
      int block_i = row / dimNvar;
      int in_i = row % dimNvar;
      
      // Fast binary-search lookup in sparse K
      double k_val = get_sparse_val(in_i, in_j, K_p_ptr, K_i_ptr, K_x_ptr);
      
      if (k_val == 0.0) continue;
      
      double val = C_x_ptr[idx] * k_val;
      
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

// C is the dense inverse of the coefficient matrix.
// K is the dense constant matrix (GRM).
// [[Rcpp::export]]
NumericMatrix getBlockTracesDenseDense(NumericVector C, int C_rows, NumericVector K, int dimVar, int dimNvar, int offset) {
  
  NumericMatrix traces(dimVar, dimVar); // Matrix of traces.
  const double* C_ptr = REAL(C); // Pointer to the data in the complete dense C^-1 (incl. fixed effects!)
  const double* K_ptr = REAL(K); // Pointer to the data in the dense K^-1.
  
  for (int i = 0; i < dimVar; ++i) { // For each row in traces or equivalently for each row in sigma^-1.
    for (int j = i; j < dimVar; ++j) { // Col: Start at the diagonal and move to the right, so only do the upper triangle of traces/sigma^-1.
      double tr = 0.0; // Initialize the trace at 0.
      
      // offset is the diagonal element (0-based indexing) of C^-1 where the block for the
      // current random effect starts. For a GBLUP model with a single fixed intercept and
      // a single random effect, offset starts at 1, C^-1[0, 0] is the intercept.
      // After processing a random effect, offset = offset + dim.var * dim.nvar.
      int rstart = offset + i * dimNvar; // Which row in C^-1 does the current block start at?
      int cstart = offset + j * dimNvar; // Which col in C^-1 does the current block start at?
      
      for (int c = 0; c < dimNvar; ++c) { // Moving along the columns of K^-1.
        // K^-1 is in col-major linear storage, so to get the index for an element in linear
        // storage, we need to add the number of rows (dimNvar) as we increase the col index:
        int k_col_offset = c * dimNvar;
        // Similarly, we need an offset for C^-1, however, here we also need to take into
        // account that we do not start at the top-left corner of C^-1. Where we start depends
        // on the number of fixed effects and what the current random effect is. For the GBLUP
        // example: cstart will be 1, meaning we already move ahead C_rows in linear storage to
        // ignore the first column of C^-1:
        int ci_col_offset = (cstart + c) * C_rows;
        
        for (int r = 0; r <= c; ++r) { // Iterate only over the upper triangle (r <= c).
          
          // This is where it all comes together to compute the exact index for the current element in
          // K^-1 and the corresponding element in C^-1. Given the col-major format, for K^-1 it's
          // simple: idx = c * nrow(K^-1) + r, using 0-based indexing and nrow(K^-1) = dimNvar.
          // For C^-1:
          // idx = rstart + r + ci_col_offset
          //     = rstart + r + (cstart + c) * C_rows
          // where rstart is the row-origin for the current sub-block of the current diagonal block of C^-1.
          // This sub-block corresponds to the product of an element of sigma^-1 with the matrix K^-1.
          // r is the row within that sub-block (rows of K^-1), cstart is the col-origin
          // for the current sub-block of the current diagonal block of C^-1, note that
          // rstart does not have to be the same as cstart, but rstart can never be
          // larger than cstart as we only do the upper triangle.
          // c is the col within the sub-block.
          double val = K_ptr[r + k_col_offset] * C_ptr[rstart + r + ci_col_offset];
          if (r == c) {
            tr += val; // Diagonal: count once
          } else {
            tr += 2.0 * val; // Off-diagonal: account for symmetry (r,c) and (c,r)
          }
        }
      }
      traces(i, j) = tr;
    }
  }
  return traces;
}

// C is the dense inverse of the coefficient matrix.
// K is the sparse constant matrix (pedigree based matrix).
// [[Rcpp::export]]
NumericMatrix getBlockTracesDenseSparse(NumericVector C, int C_rows, IntegerVector K_p, IntegerVector K_i, NumericVector K_x, int dimVar, int dimNvar, int offset) {
  
  NumericMatrix traces(dimVar, dimVar);
  const double* C_ptr = REAL(C);
  const int* K_p_ptr = INTEGER(K_p);
  const int* K_i_ptr = INTEGER(K_i);
  const double* K_x_ptr = REAL(K_x);
  
  int n_cols = K_p.length() - 1;
  for (int c = 0; c < n_cols; ++c) {
    int start = K_p_ptr[c];
    int end = K_p_ptr[c + 1];
    
    for (int idx = start; idx < end; ++idx) {
      int r = K_i_ptr[idx];
      double val = K_x_ptr[idx];
      
      for (int i = 0; i < dimVar; ++i) {
        // Global row/col indices for the stored (r, c) entry
        int row_1 = offset + i * dimNvar + r;
        
        for (int j = i; j < dimVar; ++j) {
          int col_1 = offset + j * dimNvar + c;
          
          // 1. Contribution from the stored entry (r, c)
          traces(i, j) += val * C_ptr[row_1 + col_1 * C_rows];
          
          // 2. Contribution from the symmetric counterpart (c, r) if off-diagonal
          if (r != c) {
            int row_2 = offset + i * dimNvar + c;
            int col_2 = offset + j * dimNvar + r;
            traces(i, j) += val * C_ptr[row_2 + col_2 * C_rows];
          }
        }
      }
    }
  }
  return traces;
}



// [[Rcpp::export]]
double getTraceSparseDense(IntegerVector p, IntegerVector i, NumericVector x, NumericVector dense_x, int n_rows) {
  double tr = 0.0;
  int n_cols = p.size() - 1;
  
  const int* p_ptr = INTEGER(p);
  const int* i_ptr = INTEGER(i);
  const double* x_ptr = REAL(x);
  const double* dense_x_ptr = REAL(dense_x);
  
  for (int col = 0; col < n_cols; ++col) {
    int start = p_ptr[col];
    int end = p_ptr[col + 1];
    for (int idx = start; idx < end; ++idx) {
      int row = i_ptr[idx];
      double val = x_ptr[idx] * dense_x_ptr[row + col * n_rows];
      if (row == col) {
        tr += val;
      } else {
        tr += 2.0 * val;
      }
    }
  }
  return tr;
}



