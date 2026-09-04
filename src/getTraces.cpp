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

// C is the sparse inverse coefficient matrix.
// K is the dense constant matrix (GRM) passed as a flat vector.
// [[Rcpp::export]]
NumericMatrix getBlockTracesSparseDense(IntegerVector C_p, IntegerVector C_i, NumericVector C_x, NumericVector K_x, int dimVar, int dimNvar, int offset) {
  
  NumericMatrix traces(dimVar, dimVar); // Matrix to hold the traces, size is the same as sigma.
  
  const int* C_p_ptr = INTEGER(C_p); // Column pointers of the full C^-1.
  const int* C_i_ptr = INTEGER(C_i); // Row indices for the non-zero values of C^-1, only the diagonal and upper triangle.
  const double* C_x_ptr = REAL(C_x); // Values.
  const double* K_x_ptr = REAL(K_x); // Values of the dense GRM.
  
  int block_size = dimVar * dimNvar; // Size of the diagonal block of C^-1 for the current random effect.
  int start_col = offset; // Starting col of the diagonal block of C^-1 for the current random effect (0-based).
  int end_col = offset + block_size; // End col of the diagonal block of C^-1 for the current random effect (0-based).
  
  // Only iterate through the columns of C that belong to this block:
  for (int col = start_col; col < end_col; ++col) {
    int c_rel = col - offset; // Relative column index within the diagonal block corresponding to the current random effect.
    int block_j = c_rel / dimNvar; // Integer division rounding down to determine the sub-block col index within the current diagonal block.
    int in_j = c_rel % dimNvar; // Modulo operator, remainder to determine the col index within the sub-block.
    
    // Global start and end linear indices for the non-zero values in the current column of C^-1. Note
    // that this includes only values in the upper triangle and diagonal as C is in CSC packed format.
    // It also includes non-zero values outside the diagonal block for the current random effect. That's
    // why we need the check if (row >= start_col && row < end_col) further down.
    int start = C_p_ptr[col]; 
    int end = C_p_ptr[col + 1];
    
    // The upper triangular packed CSC storage format of C^-1 ensures we never look at elements below the diagonal
    // of C^-1! So row can never be larger than col!
    for (int idx = start; idx < end; ++idx) { // Loop over all non-zero elements in this column.
      int row = C_i_ptr[idx]; // Get the global row index for the current non-zero element.
      
      // We need a check to ensure we only process elements within the diagonal block for the
      // current random effect.
      // row >= start_col ensures we only process rows below the diagonal element (start_col, start_col).
      // row < end_col ensures we only process rows above the diagonal element (end_col, end_col) which is the
      // first diagonal element of the next diagonal block (the next random effect).
      if (row >= start_col && row < end_col) {
        int r_rel = row - offset; // Relative/local row index within the current diagonal block of C^-1.
        int block_i = r_rel / dimNvar; // Integer division rounding down to get the sub-block row index within the current diagonal block.
        int in_i = r_rel % dimNvar; // Modulo operator, remainder to determine the row index within the sub-block.
        
        // LAPACK SAFETY: Always read the upper triangle of the dense K vector
        int r_k = std::min(in_i, in_j);
        int c_k = std::max(in_i, in_j);
        double k_val = K_x_ptr[r_k + c_k * dimNvar];
        double val = C_x_ptr[idx] * k_val;
        
        // Only double the off-diagonal elements if we are inside a DIAGONAL sub-block.
        // Off-diagonal sub-blocks (block_i != block_j) are fully stored in the upper 
        // triangle of C^-1, so we don't double them:
        if (block_i == block_j && in_i != in_j) {
          val *= 2.0;
        }
        traces(block_i, block_j) += val;
      }
    }
  }
  return traces;
  
  // NumericMatrix traces(dimVar, dimVar);
  // 
  // const int* C_p_ptr = INTEGER(C_p);
  // const int* C_i_ptr = INTEGER(C_i);
  // const double* C_x_ptr = REAL(C_x);
  // const double* K_x_ptr = REAL(K_x);
  // 
  // int n_cols = C_p.length() - 1;
  // for (int col = 0; col < n_cols; ++col) {
  //   int block_j = col / dimNvar;
  //   int in_j = col % dimNvar;
  //   
  //   int start = C_p_ptr[col];
  //   int end = C_p_ptr[col + 1];
  //   
  //   for (int idx = start; idx < end; ++idx) {
  //     int row = C_i_ptr[idx];
  //     int block_i = row / dimNvar;
  //     int in_i = row % dimNvar;
  //     
  //     // LAPACK SAFETY: Always read the upper triangle of the dense K vector
  //     int r_k = std::min(in_i, in_j);
  //     int c_k = std::max(in_i, in_j);
  //     
  //     double k_val = K_x_ptr[r_k + c_k * dimNvar];
  //     double val = C_x_ptr[idx] * k_val;
  //     
  //     // No * 2.0 multiplier needed; C_sub is a dgCMatrix with both triangles
  //     traces(block_j, block_i) += val;
  //   }
  // }
  // return traces;
}

// C is the sparse block of the inverse coefficient matrix.
// K is the sparse constant matrix (pedigree based matrix).
// [[Rcpp::export]]
NumericMatrix getBlockTracesSparseSparse(IntegerVector C_p, IntegerVector C_i, NumericVector C_x, IntegerVector K_p, IntegerVector K_i, NumericVector K_x, int dimVar, int dimNvar) {
  
  NumericMatrix traces(dimVar, dimVar);
  
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
      
      double k_val = get_sparse_val(in_i, in_j, K_p_ptr, K_i_ptr, K_x_ptr);
      if (k_val == 0.0) continue;
      
      double val = C_x_ptr[idx] * k_val;
      
      // No * 2.0 multiplier needed; C_sub is a dgCMatrix with both triangles
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
  const int* K_p_ptr = INTEGER(K_p); // Col pointers.
  const int* K_i_ptr = INTEGER(K_i); // Row indices.
  const double* K_x_ptr = REAL(K_x); // Values.
  
  int n_cols = K_p.length() - 1; // Number of cols in K.
  for (int c = 0; c < n_cols; ++c) { // Looping over cols in K, and equivalently the columns of the sub-block of the diagonal block of C^-1 for the current random effect.
    int start = K_p_ptr[c]; // Index for the first value in col c of K.
    int end = K_p_ptr[c + 1]; // Index for the first value in the next col of K.
    
    // Note that K is a dsCMatrix which is a CSC sparse matrix in packed form. So we
    // have a row index for each value in the upper triangle + diagonal, as well as a
    // column pointer for each column.
    for (int idx = start; idx < end; ++idx) { // Looping over all values in col c of K.
      int r = K_i_ptr[idx]; // The row index for the current value within K, and equivalently within each sub-block of the diagonal block of C^-1 for the current random effect.
      double val = K_x_ptr[idx]; // Getting the actual value from K.
      
      for (int i = 0; i < dimVar; ++i) { // Looping over the rows of sigma / sub-blocks of the current diagonal block of C^-1.
        // Global row index within the full C^-1.
        // offset determines the starting row for the current diagonal block of C^-1.
        // i * dimNvar moves the index to the starting position of the ith sub-block of the
        // current diagonal block of C^-1.
        // r moves the index to the current row within that sub-block.
        int row_1 = offset + i * dimNvar + r; 
        
        // Looping over the cols of sigma / sub-blocks of the current diagonal block of C^-1.
        // Note that we let i define the sub-block where we start. For the DenseDense function
        // we go from top to the diagonal, then left to right.
        // For this DenseSparse function we are going from the diagonal to the right within a row of sigma,
        // then step down and to the right to the next diagonal sub-block. That's fine.
        for (int j = i; j < dimVar; ++j) {
          
          // Global col index within the full C^-1.
          // offset determines the starting col for the current diagonal block of C^-1.
          // j * dimNvar then moves the index to the starting column of the current sub-block
          // within the current diagonal block of C^-1.
          // c then moves the index to the current column within that sub-block.
          int col_1 = offset + j * dimNvar + c;
          
          // 1. Contribution from the stored entry (r, c)
          // Here the linear index for the value of C^-1 corresponding to the current
          // non-zero element of K^-1 is computed as:
          // row_1 + col_1 * C_rows.
          // C^-1 is col-major and dense, so it's a simple formula.
          if (r != c) {
            traces(i, j) += 2 * val * C_ptr[row_1 + col_1 * C_rows];
          } else {
            traces(i, j) += val * C_ptr[row_1 + col_1 * C_rows];
          }
          // 2. Contribution from the symmetric counterpart (c, r) if off-diagonal
          // if (r != c) {
          //   int row_2 = offset + i * dimNvar + c;
          //   int col_2 = offset + j * dimNvar + r;
          //   traces(i, j) += val * C_ptr[row_2 + col_2 * C_rows];
          // }
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



