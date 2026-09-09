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
// K is the dense constant matrix (inverse GRM) passed as a flat vector.
// [[Rcpp::export]]
NumericMatrix getBlockTracesSparseDense(IntegerVector C_p, IntegerVector C_i, NumericVector C_x,
                                        NumericVector K_x,
                                        int dimVar, int dimNvar, int offset) {
  
  NumericMatrix traces(dimVar, dimVar); // Matrix to hold the traces, size is the same as sigma.
  
  const int* C_p_ptr = INTEGER(C_p); // Column pointers of the full C^-1.
  const int* C_i_ptr = INTEGER(C_i); // Row indices for the non-zero values of C^-1, only the diagonal and upper triangle.
  const double* C_x_ptr = REAL(C_x); // Values.
  const double* K_x_ptr = REAL(K_x); // Values of the dense inverse GRM.
  
  int block_size = dimVar * dimNvar; // Size of the diagonal block of C^-1 for the current random effect.
  int start_col = offset; // Starting col of the diagonal block of C^-1 for the current random effect (0-based).
  int end_col = offset + block_size; // End col of the diagonal block of C^-1 for the current random effect (0-based).
  
  // Only iterate through the columns of C that belong to this block:
  for (int col = start_col; col < end_col; ++col) {
    int c_rel = col - offset; // Relative column index within the diagonal block corresponding to the current random effect.
    int block_j = c_rel / dimNvar; // Integer division rounding down to determine the sub-block col index within the current diagonal block.
    int in_j = c_rel % dimNvar; // Modulo operator, remainder to determine the col index within the sub-block.
    
    // Global start and end linear indices for the non-zero values in the current column of C^-1. Note
    // that this includes only values in the LOWER triangle and diagonal as C is in CSC packed format.
    // SIS() returns a packed CSC matrix in LOWER format!
    // It also includes non-zero values outside the diagonal block for the current random effect. That's
    // why we need the check if (row >= start_col && row < end_col) further down.
    int start = C_p_ptr[col]; 
    int end = C_p_ptr[col + 1];
    
    // The lower triangular packed CSC storage format of C^-1 ensures we never look at elements above the diagonal
    // of C^-1! So row can never be smaller than col!
    for (int idx = start; idx < end; ++idx) { // Loop over all non-zero elements in this column.
      int row = C_i_ptr[idx]; // Get the global row index for the current non-zero element.
      
      // We need a check to ensure we only process elements within the diagonal block for the
      // current random effect.
      // row >= start_col ensures we only process rows at or below the diagonal element (start_col, start_col).
      // row < end_col ensures we only process rows above the diagonal element (end_col, end_col) which is the
      // first diagonal element of the next diagonal block (the next random effect).
      if (row >= start_col && row < end_col) {
        int r_rel = row - offset; // Relative/local row index within the current diagonal block of C^-1.
        int block_i = r_rel / dimNvar; // Integer division rounding down to get the sub-block row index within the current diagonal block.
        int in_i = r_rel % dimNvar; // Modulo operator, remainder to determine the row index within the sub-block.
        
        // Always read the upper triangle of K, even if we are in an off-diagonal block so that
        // in_i can be larger than in_j:
        int r_k = std::min(in_i, in_j);
        int c_k = std::max(in_i, in_j);
        double k_val = K_x_ptr[r_k + c_k * dimNvar];
        double val = C_x_ptr[idx] * k_val;
        
        // Only double the off-diagonal elements if we are inside a DIAGONAL sub-block.
        // Off-diagonal sub-blocks (block_i != block_j) are fully stored in the lower 
        // triangle of C^-1, so we don't double them:
        if (block_i == block_j && in_i != in_j) {
          val *= 2.0;
        }
        // If we put traces(block_i, block_j) here we would be writing into the
        // lower triangle of traces because the SIS is in lower packed format.
        // I want to write into the upper triangle to keep it consistent with the
        // DenseDense and DenseSparse functions, so we simply swap i and j here:
        traces(block_j, block_i) += val;
      }
    }
  }
  return traces;
}

// C is the sparse inverse coefficient matrix.
// K is the sparse constant matrix (inverse pedigree based matrix).
// [[Rcpp::export]]
NumericMatrix getBlockTracesSparseSparse(IntegerVector C_p, IntegerVector C_i, NumericVector C_x, 
                                         IntegerVector K_p, IntegerVector K_i, NumericVector K_x, 
                                         int dimVar, int dimNvar, int offset) {
  
  NumericMatrix traces(dimVar, dimVar); // Output matrix of same dimensions as sigma.
  
  const int* C_p_ptr = INTEGER(C_p);
  const int* C_i_ptr = INTEGER(C_i);
  const double* C_x_ptr = REAL(C_x);
  const int* K_p_ptr = INTEGER(K_p);
  const int* K_i_ptr = INTEGER(K_i);
  const double* K_x_ptr = REAL(K_x);
  
  int block_size = dimVar * dimNvar; // Size of the diagonal block of C^-1 for the current random effect.
  int start_col = offset; // Where does the current diagonal block of C^-1 start?
  int end_col = offset + block_size; // Where does it end?
  
  // Only iterate through the cols of C^-1 that belong to the current diagonal block:
  for (int col = start_col; col < end_col; ++col) {
    int c_rel = col - offset; // Column within the current diagonal block.
    int block_j = c_rel / dimNvar; // Col-coordinate of the current sub-block.
    int in_j = c_rel % dimNvar; // Local column number within the sub-block.
    
    int start = C_p_ptr[col]; // Start index of the values in this col of C^-1.
    int end = C_p_ptr[col + 1]; // Start index of the values in the next col of C^-1.
    
    // C^-1 is lower triangular packed CSC from SIS(), so row >= col is guaranteed:
    for (int idx = start; idx < end; ++idx) { // Iterating over the non-zero values in the current col of C^-1.
      int row = C_i_ptr[idx]; // Global row index of that value within C^-1.
      
      // Check so we only process values within the current diagonal block of C^-1.
      // If row < start_col we are above the diagonal block, if row >= end_col we
      // are below the diagonal block:
      if (row >= start_col && row < end_col) {
        int r_rel = row - offset; // Row within the current diagonal block of C^-1.
        int block_i = r_rel / dimNvar; // Row coordinate of the current sub-block.
        int in_i = r_rel % dimNvar; // Local row number within the sub-block.
        
        // Fetch k_val using the helper. It is 100% agnostic to K's uplo and 
        // perfectly handles lookups where in_i > in_j in asymmetric off-diagonal blocks.
        double k_val = get_sparse_val(in_i, in_j, K_p_ptr, K_i_ptr, K_x_ptr);
        if (k_val == 0.0) { // No need to do anything.
          continue;
        }
        double val = C_x_ptr[idx] * k_val;
        
        // C^-1 logic is identical to SparseDense.
        // Only double the off-diagonal elements if we are inside a diagonal sub-block.
        // Off-diagonal sub-blocks (block_i != block_j) are fully stored in the lower 
        // triangle of C^-1, so we don't double them:
        if (block_i == block_j && in_i != in_j) {
          val *= 2.0;
        }
        // Swap indices to write into the upper triangle of traces:
        traces(block_j, block_i) += val;
      }
    }
  }
  return traces;
}

// C is the dense inverse of the coefficient matrix.
// K is the dense constant matrix (inverse GRM).
// [[Rcpp::export]]
NumericMatrix getBlockTracesDenseDense(NumericVector C, int C_rows,
                                       NumericVector K,
                                       int dimVar, int dimNvar, int offset) {
  
  NumericMatrix traces(dimVar, dimVar); // Output matrix, same size as sigma for the current random effect.
  const double* C_ptr = REAL(C); // Values in dense C^-1.
  const double* K_ptr = REAL(K); // Values in dense K^-1.
  
  for (int i = 0; i < dimVar; ++i) { // Iterating over the rows of traces/sigma.
    for (int j = i; j < dimVar; ++j) { // Iterating over the diagonal + cols to the right for the current row of traces/sigma
      double tr = 0.0;
      // Below, offset moves the coordinate to the diagonal block for the current random effect.
      // i/j * dimNvar then moves the coordinate to the current subblock corresponding to an element
      // of traces/sigma:
      int rstart = offset + i * dimNvar; // Global index for the starting row.
      int cstart = offset + j * dimNvar; // Global index for the starting col.
      
      // Iterating over the columns of K^-1:
      for (int c = 0; c < dimNvar; ++c) {
        // Same thing for C^-1, but we also need to keep the offset and sub-blocks in mind:
        int ci_col_offset = (cstart + c) * C_rows;
        
        if (i == j) {
          // Diagonal sub-block: always symmetric and only upper half-stored.
          // Only read the upper triangle (r <= c):
          for (int r = 0; r <= c; ++r) { // For a given col, start at the first row and move down to the diagonal.
            double val = K_ptr[r + c * dimNvar] * C_ptr[rstart + r + ci_col_offset]; // Compute value.
            if (r == c) {
              tr += val; // Add diagonal values once.
            } else {
              tr += 2.0 * val; // Add off-diagonal values twice to account for symmetry.
            }
          }
        } else {
          // Off-diagonal sub-block: not always symmetric, so we need to read the full sub-block.
          // This is fine, because the off-diagonal sub-block (i, j) exists fully in the upper triangle
          // of C^-1. So it is fully stored.
          for (int r = 0; r < dimNvar; ++r) { // For a given col, start at the first row and move all the way down.
            // K^-1 is a dsyMatrix, which I think only stores meaningful values in the upper triangle.
            // The lower triangle does exist in memory (dsyMatrix is not packed), but it might contain
            // nonsense. Just to be safe, only read from the upper triangle of K^-1.
            int r_k = std::min(r, c);
            int c_k = std::max(r, c);
            double k_val = K_ptr[r_k + c_k * dimNvar];
            // No doubling needed because (1) the sub-block is not guaranteed to be symmetric,
            // and (2) we read the full thing anyway:
            tr += k_val * C_ptr[rstart + r + ci_col_offset];
          }
        }
      }
      traces(i, j) = tr;
    }
  }
  return traces;
}

// C is the dense inverse of the coefficient matrix.
// K is the sparse constant matrix (inverse pedigree based matrix).
NumericMatrix getBlockTracesDenseSparse(NumericVector C, int C_rows,
                                        IntegerVector K_p, IntegerVector K_i, NumericVector K_x,
                                        int dimVar, int dimNvar, int offset) {

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
        for (int j = i; j < dimVar; ++j) {

          // Global col index within the full C^-1.
          // offset determines the starting col for the current diagonal block of C^-1.
          // j * dimNvar then moves the index to the starting column of the current sub-block
          // within the current diagonal block of C^-1.
          // c then moves the index to the current column within that sub-block.
          int col_1 = offset + j * dimNvar + c;
          
          if (i == j) { // If we are dealing with a diagonal sub-block.
            if (r != c) { // If we are in the upper triangle of K^-1 (off-diagonal)
              traces(i, j) += 2 * val * C_ptr[row_1 + col_1 * C_rows];
            } else {
              traces(i, j) += val * C_ptr[row_1 + col_1 * C_rows];
            }
          } else { // If i < j, so an off-diagonal sub-block in the upper triangle of C^-1.
            if (r == c) {
              traces(i, j) += val * C_ptr[row_1 + col_1 * C_rows];
            } else { // Do contributions from the upper and lower triangle separately:
              // Upper triangle of the off-diagonal sub-block of C^-1:
              traces(i, j) += val * C_ptr[row_1 + col_1 * C_rows];
              // Lower triangle of the off-diagonal sub-block of C^-1 (not necessarily the
              // same as the value in the upper triangle!):
              int row_2 = offset + i * dimNvar + c;
              int col_2 = offset + j * dimNvar + r;
              traces(i, j) += val * C_ptr[row_2 + col_2 * C_rows];
            }
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



