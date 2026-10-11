#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// In-place update of the sparse inverse subset:
//
//   Ciss[i, j] += (Q %*% Cidd %*% t(Q))[i, j]   for every *stored* entry (i, j) of Ciss
//
// Ciss  : dsCMatrix (M x M), either triangle. Modified in place.
// QCidd : dgeMatrix (M x N), column-major, QCidd = Q %*% Cidd. Computed on the CPU (Q %*% Cidd) or on the GPU
//         (cuSolverPhase2(..., Qt)$QCidd); the layout is identical.
//
// Both functions below are serial and compute the same thing; they differ only in loop order (and in whether they
// take Q or t(Q)). Per entry they add the terms in increasing k, so they give bit-identical results.
//
//   Ciss[i, j] += sum_k QCidd[i, k] * Q[j, k]


// ---- Variant 1 (k outermost, takes Q) ---------------------------------------------------------------------------
// Q Cidd Q' = sum_k QCidd[, k] %o% Q[, k]: an outer product of a *contiguous* column of QCidd with a *sparse* column
// of Q (a CSC column of Q), sampled at the stored entries of Ciss. With k outermost, column k of QCidd (M doubles,
// ~170 KB for M = 21306) stays hot in cache for the whole pass over the nonzeros of Q[, k], and the arrays of Ciss
// (x and i, ~6.6 MB here) are small enough to live in the outer cache levels. The big M x N array is only ever read
// one column at a time, front to back.
//
// Q : dgCMatrix (M x N)
// [[Rcpp::export]]
void inplaceUpdateSIS(S4 Ciss, S4 QCidd, S4 Q) {
  NumericVector cx   = Ciss.slot("x");      // modified in place, no copy
  IntegerVector cp   = Ciss.slot("p");
  IntegerVector ci   = Ciss.slot("i");
  IntegerVector cdim = Ciss.slot("Dim");
  
  NumericVector X    = QCidd.slot("x");
  IntegerVector xdim = QCidd.slot("Dim");
  
  NumericVector qx   = Q.slot("x");
  IntegerVector qp   = Q.slot("p");
  IntegerVector qi   = Q.slot("i");
  IntegerVector qdim = Q.slot("Dim");
  
  const int M = cdim[0];
  const int N = qdim[1];
  if (cdim[1] != M || qdim[0] != M || xdim[0] != M || xdim[1] != N) {
    stop("Dimension mismatch: Ciss is %i x %i, QCidd is %i x %i, Q is %i x %i (expected M x M, M x N, M x N).",
         cdim[0], cdim[1], xdim[0], xdim[1], qdim[0], qdim[1]);
  }
  
  double*       cx_ptr = cx.begin();
  const int*    cp_ptr = cp.begin();
  const int*    ci_ptr = ci.begin();
  const double* X_ptr  = X.begin();
  const double* qx_ptr = qx.begin();
  const int*    qp_ptr = qp.begin();
  const int*    qi_ptr = qi.begin();
  
  for (int k = 0; k < N; ++k) {
    const double* Xk = X_ptr + static_cast<R_xlen_t>(k) * M;   // contiguous column k of QCidd
    
    for (int t = qp_ptr[k]; t < qp_ptr[k + 1]; ++t) {
      const int    j     = qi_ptr[t];                          // row of Q = column of Ciss
      const double w     = qx_ptr[t];                          // Q[j, k]
      const int    c_end = cp_ptr[j + 1];
      
      for (int idx = cp_ptr[j]; idx < c_end; ++idx) {
        cx_ptr[idx] += w * Xk[ci_ptr[idx]];                    // QCidd[i, k] * Q[j, k]
      }
    }
  }
}


// ---- Variant 2 (j outermost, takes Qt) --------------------------------------------------------------------------
// Column j of Qt is row j of Q, so Ciss[, j] only needs the nonzeros of Qt[, j]. Each (j, k) pair then reads a
// different contiguous column k of QCidd at the row indices of Ciss[, j]. Kept mainly as a baseline to benchmark
// against variant 1 (and it is the form that parallelises over j without races if you ever want OpenMP).
//
// Qt : dgCMatrix (N x M), Qt = t(Q)
// [[Rcpp::export]]
void inplaceUpdateSIS_Qt(S4 Ciss, S4 QCidd, S4 Qt) {
  NumericVector cx   = Ciss.slot("x");
  IntegerVector cp   = Ciss.slot("p");
  IntegerVector ci   = Ciss.slot("i");
  IntegerVector cdim = Ciss.slot("Dim");
  
  NumericVector X    = QCidd.slot("x");
  IntegerVector xdim = QCidd.slot("Dim");
  
  NumericVector tx   = Qt.slot("x");
  IntegerVector tp   = Qt.slot("p");
  IntegerVector ti   = Qt.slot("i");
  IntegerVector tdim = Qt.slot("Dim");
  
  const int M = cdim[0];
  const int N = tdim[0];
  if (cdim[1] != M || tdim[1] != M || xdim[0] != M || xdim[1] != N) {
    stop("Dimension mismatch: Ciss is %i x %i, QCidd is %i x %i, Qt is %i x %i (expected M x M, M x N, N x M).",
         cdim[0], cdim[1], xdim[0], xdim[1], tdim[0], tdim[1]);
  }
  
  double*       cx_ptr = cx.begin();
  const int*    cp_ptr = cp.begin();
  const int*    ci_ptr = ci.begin();
  const double* X_ptr  = X.begin();
  const double* tx_ptr = tx.begin();
  const int*    tp_ptr = tp.begin();
  const int*    ti_ptr = ti.begin();
  
  for (int j = 0; j < M; ++j) {
    const int c_start = cp_ptr[j];
    const int c_end   = cp_ptr[j + 1];
    if (c_start == c_end) continue;
    
    for (int t = tp_ptr[j]; t < tp_ptr[j + 1]; ++t) {
      const int     k  = ti_ptr[t];                              // column of QCidd
      const double  w  = tx_ptr[t];                              // Q[j, k]
      const double* Xk = X_ptr + static_cast<R_xlen_t>(k) * M;
      
      for (int idx = c_start; idx < c_end; ++idx) {
        cx_ptr[idx] += w * Xk[ci_ptr[idx]];
      }
    }
  }
}


// [[Rcpp::export]]
void inplaceUpdateSIS_old(S4 Ciss, S4 tW_mat, S4 tQ_mat) {
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