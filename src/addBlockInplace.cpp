#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void addBlockInplace(S4 C, S4 Gi, int offset) {
  // 1. Extract the SEXP objects
  NumericVector Cx = C.slot("x");
  NumericVector Gix = Gi.slot("x");
  
  // 2. Extract RAW C++ double pointers. 
  // This entirely bypasses Rcpp's [] operator overhead.
  double* ptr_C = REAL(Cx);
  double* ptr_Gi = REAL(Gix);
  
  IntegerVector C_dim = C.slot("Dim");
  IntegerVector Gi_dim = Gi.slot("Dim");
  
  // 3. Use size_t for fast memory addressing and to prevent 32-bit overflow
  size_t N_C = C_dim[0];
  size_t N_Gi = Gi_dim[0];
  size_t off = offset - 1; 
  
  for (size_t j = 0; j < N_Gi; j++) {
    // 4. Hoist the heavy math outside the inner loop!
    // Calculate the base index for the current column just ONCE per column.
    size_t col_C = (j + off) * N_C + off;
    size_t col_Gi = j * N_Gi;
    
    for (size_t i = 0; i < N_Gi; i++) {
      // 5. Bare-metal pointer arithmetic
      ptr_C[col_C + i] += ptr_Gi[col_Gi + i];
    }
  }
}