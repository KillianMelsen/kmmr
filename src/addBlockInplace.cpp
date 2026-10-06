#include <Rcpp.h>
#include <algorithm> // Required for std::fill
using namespace Rcpp;

// [[Rcpp::export]]
void resetCoefficientMatrix(S4 base) {
  // Extract slots and reset numerical values to 0.0 while preserving the sparse pattern:
  NumericVector base_x_ptr = base.slot("x");
  std::fill(base_x_ptr.begin(), base_x_ptr.end(), 0.0);
  
  // Clear the cached Cholesky factor! Otherwise any Cholesky() call in R will just
  // return these factors which are now outdated since we have updated the values
  // in-place (which R does not know about):
  if (base.hasSlot("factors")) {
    base.slot("factors") = Rcpp::List::create();
  }
}


// [[Rcpp::export]]
void addToCoefficientMatrix(S4 base, S4 addition, int offset) {
  
  bool base_sparse = base.hasSlot("p");
  bool add_sparse = addition.hasSlot("p");
  
  if (base_sparse && add_sparse) { 
    
    // ==========================================
    // 1. SPARSE BASE + SPARSE ADDITION
    // ==========================================
    IntegerVector base_p = base.slot("p");
    IntegerVector base_i = base.slot("i");
    NumericVector base_x = base.slot("x");
    
    IntegerVector add_p = addition.slot("p");
    IntegerVector add_i = addition.slot("i");
    NumericVector add_x = addition.slot("x");
    
    int* base_p_ptr = INTEGER(base_p);
    int* base_i_ptr = INTEGER(base_i);
    double* base_x_ptr = REAL(base_x);
    
    int* add_p_ptr = INTEGER(add_p);
    int* add_i_ptr = INTEGER(add_i);
    double* add_x_ptr = REAL(add_x);
    
    int add_cols = add_p.size() - 1;
    
    for (int j = 0; j < add_cols; ++j) {
      int add_start = add_p_ptr[j];
      int add_end   = add_p_ptr[j + 1];
      if (add_start == add_end) { continue; }
      
      int base_col = j + offset;
      int base_idx = base_p_ptr[base_col];
      int base_end = base_p_ptr[base_col + 1];
      
      for (int k = add_start; k < add_end; ++k) {
        int target_row = add_i_ptr[k] + offset;
        
        while (base_idx < base_end && base_i_ptr[base_idx] < target_row) {
          base_idx++;
        }
        
        if (base_idx < base_end && base_i_ptr[base_idx] == target_row) {
          base_x_ptr[base_idx] += add_x_ptr[k]; 
        } else {
          Rcpp::stop("Sparsity mismatch: Element in addition matrix not found in base matrix!");
        }
      }
    }
    
  } else if (!base_sparse && !add_sparse) { 
    
    // ==========================================
    // 2. DENSE BASE + DENSE ADDITION
    // ==========================================
    NumericVector base_x = base.slot("x");
    NumericVector add_x = addition.slot("x");
    
    double* base_x_ptr = REAL(base_x);
    double* add_x_ptr = REAL(add_x);
    
    IntegerVector base_dim = base.slot("Dim");
    IntegerVector add_dim = addition.slot("Dim");
    
    size_t base_n = base_dim[0];
    size_t add_n = add_dim[0];
    size_t offset64 = offset; 
    
    for (size_t j = 0; j < add_n; j++) {
      size_t add_col = j * add_n; 
      size_t base_col = (j + offset64) * base_n + offset64; 
      
      for (size_t i = 0; i < add_n; i++) { 
        base_x_ptr[base_col + i] += add_x_ptr[add_col + i];
      }
    }
    
  } else if (!base_sparse && add_sparse) {
    
    // ==========================================
    // 3. DENSE BASE + SPARSE ADDITION 
    // ==========================================
    NumericVector base_x = base.slot("x");
    double* base_x_ptr = REAL(base_x);
    
    IntegerVector base_dim = base.slot("Dim");
    size_t base_n = base_dim[0];
    size_t offset64 = offset; 
    
    IntegerVector add_p = addition.slot("p");
    IntegerVector add_i = addition.slot("i");
    NumericVector add_x = addition.slot("x");
    
    int* add_p_ptr = INTEGER(add_p);
    int* add_i_ptr = INTEGER(add_i);
    double* add_x_ptr = REAL(add_x);
    
    int add_cols = add_p.size() - 1;
    
    for (int j = 0; j < add_cols; ++j) {
      int add_start = add_p_ptr[j];
      int add_end   = add_p_ptr[j + 1];
      
      // Calculate the memory start of this column in the dense base matrix
      size_t base_col_start = (j + offset64) * base_n;
      
      for (int k = add_start; k < add_end; ++k) {
        size_t target_row = add_i_ptr[k] + offset64;
        
        // Add the sparse element exactly where it belongs in the dense grid
        base_x_ptr[base_col_start + target_row] += add_x_ptr[k];
      }
    }
    
  } else {
    // 4. SPARSE BASE + DENSE ADDITION (Safety Net)
    Rcpp::stop("Cannot dynamically add a dense matrix into a sparse base matrix in-place.");
  }
}








