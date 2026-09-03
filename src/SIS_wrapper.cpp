#include <Rcpp.h>
#include <Matrix.h>
#include <Matrix/cholmod-utils.h>
#include <Matrix/alloca.h>
#include "cholmod_extra.h"
#include "SIS_impl.h"

// [[Rcpp::depends(Matrix)]]
using namespace Rcpp;

//' Compute a Sparse Inverse Subset (SIS) from a sparse Cholesky factor (CF).
//'
//' @param choleskyFactor A sparse Cholesky factor (dCHMsimpl or dCHMsuper) from Matrix::Cholesky()
// [[Rcpp::export]]
SEXP SIS(SEXP choleskyFactor) {
 
 // Initialize the CHOLMOD common environment:
 cholmod_common c;
 M_cholmod_start(&c);
 
 // Convert R object (dCHMsuper or dCHMsimpl) to a cholmod_factor pointer:
 CHM_FR L = AS_CHM_FR(choleskyFactor);
 
 // Compute the sparse inverse subset:
 cholmod_sparse *SIS = NULL;
 SIS = cholmod_SIS(L, &c);
 
 if (SIS == NULL) {
   M_cholmod_finish(&c);
   stop("cholmod_SIS failed to compute the sparse inverse subset.");
 }
 
 // Convert the resulting cholmod_sparse back to an R dsCMatrix:
 SEXP res = M_chm_sparse_to_SEXP(SIS, 1, 0, 0, "", R_NilValue);
 
 // Clean up the CHOLMOD environment:
 M_cholmod_finish(&c);
 
 return res;
}
 
 
 