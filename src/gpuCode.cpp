// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cusparse.h> // NEW: remember to link with -lcusparse (PKG_LIBS in Makevars)

#define CHECK_CUDA(call) {                                                                          \
cudaError_t err = call;                                                                             \
if (err != cudaSuccess) {                                                                           \
  Rcpp::stop("CUDA error at line %d of /src/gpuCode.cpp: %s", __LINE__, cudaGetErrorString(err));   \
}                                                                                                   \
}

// NEW:
#define CHECK_CUSPARSE(call) {                                                                                \
cusparseStatus_t st_ = call;                                                                                  \
if (st_ != CUSPARSE_STATUS_SUCCESS) {                                                                         \
  Rcpp::stop("cuSPARSE error at line %d of /src/gpuCode.cpp: %s", __LINE__, cusparseGetErrorString(st_));     \
}                                                                                                             \
}

// NEW:
#define CHECK_CUBLAS(call) {                                                                                  \
cublasStatus_t st_ = call;                                                                                    \
if (st_ != CUBLAS_STATUS_SUCCESS) {                                                                           \
  Rcpp::stop("cuBLAS error at line %d of /src/gpuCode.cpp: status code %d", __LINE__, static_cast<int>(st_)); \
}                                                                                                             \
}

// NEW: Grow-only device buffer. Only used for the cuSPARSE workspace, whose size comes from cusparseSpMM_bufferSize.
struct DevBuf {
  void*  ptr = nullptr;
  size_t cap = 0; // capacity in bytes
  
  void reserve(size_t bytes) {
    if (bytes <= cap) return;
    release();
    CHECK_CUDA( cudaMalloc(&ptr, bytes) );
    cap = bytes;
  }
  
  void release() {
    if (ptr) cudaFree(ptr);
    ptr = nullptr;
    cap = 0;
  }
};

// The Persistent Context Structure:
struct cuSolverContext {
  cublasHandle_t cublasH;
  cusolverDnHandle_t cusolverH;
  cusparseHandle_t cusparseH; // NEW
  cudaStream_t stream;
  
  int N;
  int rhs_cols_1; // Number of columns of the phase 1 right-hand side, always 1
  int rhs_cols_2; // Number of columns of the phase 2 right-hand side, equal to the number of parameters
  int lwork; // Stores the max workspace size needed
  
  // Persistent memory blocks
  double *d_C;
  double *d_RHS;
  double *d_buWu;
  double *d_work;
  double *d_diag;
  int *d_info;
  
  // NEW: everything for QCidd = Q %*% Cidd. M == 0 means Q is not used.
  // M is known at cuSolverInit, so d_Qp, d_QCidd, descB and descC are created there. The number of stored entries
  // of Q is only known at the first phase 2 call with Q; from then on it is fixed (nnzQ != 0) and d_Qi, d_Qx and
  // descA are created at that point.
  int M;                      // rows of Q (= dimension of Ciss)
  int nnzQ;                   // stored entries of Q; 0 until the first phase 2 call with Q, fixed afterwards
  int *d_Qp;                  // row pointers of Q (= column pointers of Qt), M + 1 ints
  int *d_Qi;                  // column indices of Q (= row indices of Qt), nnzQ ints
  double *d_Qx;               // values of Q, nnzQ doubles
  double *d_QCidd;            // dense M x N result, column-major
  cusparseSpMatDescr_t descA; // Q: M x N CSR view of d_Qp / d_Qi / d_Qx
  cusparseDnMatDescr_t descB; // Cidd: N x N view of d_C
  cusparseDnMatDescr_t descC; // QCidd: M x N view of d_QCidd
  DevBuf d_spmmBuf;           // cuSPARSE workspace (grow-only)
};

// The finalizer (safely tears down the struct and VRAM):
void context_finalizer_hook(SEXP ptr_in) {
  if (TYPEOF(ptr_in) == EXTPTRSXP) {
    cuSolverContext* ctx = static_cast<cuSolverContext*>(R_ExternalPtrAddr(ptr_in));
    if (ctx != nullptr) {
      if (ctx->d_C) cudaFree(ctx->d_C);
      if (ctx->d_RHS) cudaFree(ctx->d_RHS);
      if (ctx->d_buWu) cudaFree(ctx->d_buWu);
      if (ctx->d_work) cudaFree(ctx->d_work);
      if (ctx->d_diag) cudaFree(ctx->d_diag);
      if (ctx->d_info) cudaFree(ctx->d_info);
      if (ctx->d_Qp) cudaFree(ctx->d_Qp);                 // NEW
      if (ctx->d_Qi) cudaFree(ctx->d_Qi);                 // NEW
      if (ctx->d_Qx) cudaFree(ctx->d_Qx);                 // NEW
      if (ctx->d_QCidd) cudaFree(ctx->d_QCidd);           // NEW
      ctx->d_spmmBuf.release();                           // NEW
      if (ctx->descA) cusparseDestroySpMat(ctx->descA);   // NEW
      if (ctx->descB) cusparseDestroyDnMat(ctx->descB);   // NEW
      if (ctx->descC) cusparseDestroyDnMat(ctx->descC);   // NEW
      if (ctx->cublasH) cublasDestroy(ctx->cublasH);
      if (ctx->cusolverH) cusolverDnDestroy(ctx->cusolverH);
      if (ctx->cusparseH) cusparseDestroy(ctx->cusparseH); // NEW
      if (ctx->stream) cudaStreamDestroy(ctx->stream);
      
      delete ctx; // Free the struct itself
      R_ClearExternalPtr(ptr_in); // Prevent double-free
    }
  }
}

// The initializer (called from R once before iterating):
// CHANGED: optional M (= nrow(Q)) lets us allocate the dense QCidd buffer and its descriptors up front.
// Leave it at 0 if Q is not used.
// [[Rcpp::export(cuSolverInit)]]
SEXP cuSolverInit(int N, int rhs_cols_1, int rhs_cols_2, int M = 0) {
  cuSolverContext* ctx = new cuSolverContext(); // value-initialised: all raw pointers/handles start as nullptr
  ctx->N = N;
  ctx->rhs_cols_1 = rhs_cols_1; // always 1
  ctx->rhs_cols_2 = rhs_cols_2; // number of parameters (variance components)
  
  cublasCreate(&(ctx->cublasH));
  cusolverDnCreate(&(ctx->cusolverH));
  CHECK_CUDA( cudaStreamCreateWithFlags(&(ctx->stream), cudaStreamNonBlocking) );
  cublasSetStream(ctx->cublasH, ctx->stream);
  cusolverDnSetStream(ctx->cusolverH, ctx->stream);
  CHECK_CUSPARSE( cusparseCreate(&(ctx->cusparseH)) );          // NEW
  CHECK_CUSPARSE( cusparseSetStream(ctx->cusparseH, ctx->stream) ); // NEW
  
  // Allocate all memory upfront:
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_C)), N * N * sizeof(double)) );
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_RHS)), N * rhs_cols_1 * sizeof(double)) );
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_buWu)), N * rhs_cols_2 * sizeof(double)) );
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_diag)), N * sizeof(double)) );
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_info)), sizeof(int)) );
  
  // Calculate max lwork needed for both Phase 1 (potrf) and Phase 2 (potri)
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
  int lwork_potrf = 0;
  int lwork_potri = 0;
  cusolverDnDpotrf_bufferSize(ctx->cusolverH, uplo, N, ctx->d_C, N, &lwork_potrf);
  cusolverDnDpotri_bufferSize(ctx->cusolverH, uplo, N, ctx->d_C, N, &lwork_potri);
  
  ctx->lwork = std::max(lwork_potrf, lwork_potri);
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_work)), ctx->lwork * sizeof(double)) );
  
  // NEW: what we can already allocate for QCidd = Q %*% Cidd. B (Cidd) lives in d_C and C (QCidd) in d_QCidd; both
  // pointers never change, so their cuSPARSE descriptors are created once here and reused by every phase 2.
  if (M > 0) {
    ctx->M = M;
    CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_Qp)), (static_cast<size_t>(M) + 1) * sizeof(int)) );
    CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_QCidd)), static_cast<size_t>(M) * N * sizeof(double)) );
    CHECK_CUSPARSE( cusparseCreateDnMat(&(ctx->descB), N, N, N, ctx->d_C, CUDA_R_64F, CUSPARSE_ORDER_COL) );
    CHECK_CUSPARSE( cusparseCreateDnMat(&(ctx->descC), M, N, M, ctx->d_QCidd, CUDA_R_64F, CUSPARSE_ORDER_COL) );
  }
  
  // Wrap and return to R
  Rcpp::XPtr<cuSolverContext> ptr_out(ctx, false);
  R_RegisterCFinalizerEx(ptr_out, context_finalizer_hook, TRUE);
  
  return ptr_out;
}

// Manual cleanup hook for R
// [[Rcpp::export(cuSolverCleanup)]]
void cuSolverCleanup(SEXP ptr_in) {
  context_finalizer_hook(ptr_in);
}

// Device properties function
//[[Rcpp::export(.getDeviceProps)]]
Rcpp::List getDeviceProps(const int& id) {
  
  // Create CUDA device properties variable
  cudaDeviceProp prop;
  
  // Query GPU properties
  CHECK_CUDA( cudaGetDeviceProperties(&prop, id) );
  
  // Storing and returning results
  char deviceName[256];
  strcpy(deviceName, prop.name);
  return Rcpp::List::create(Rcpp::Named("id") = id,
                            Rcpp::Named("deviceName") = deviceName,
                            Rcpp::Named("integr") = prop.integrated,
                            Rcpp::Named("mjr") = prop.major,
                            Rcpp::Named("mnr") = prop.minor);
}

// Helper function to extract and cast to dmat
arma::dmat convert_to_dmat(SEXP x) {
  if (Rf_isS4(x)) {
    Rcpp::S4 obj(x);
    if (obj.hasSlot("x") && obj.hasSlot("Dim")) {
      Rcpp::NumericVector data = obj.slot("x");
      Rcpp::IntegerVector dims = obj.slot("Dim");
      
      // Initialize dmat instantly using the memory pointer:
      return arma::dmat(data.begin(), dims[0], dims[1], false, true);
    }
  }
  Rcpp::stop("Input must be a dense Matrix object (dgeMatrix or dsyMatrix)...");
}

// Mirrors the upper triangle of the column-major n x n matrix d_A (on the device) into its lower triangle.
// After potri only the upper triangle of d_C holds C^-1 (the lower triangle is leftover input), but cuSPARSE has no
// notion of a symmetric / uplo dense operand, so SpMM needs the full matrix.
// Element (r, c) lives at d_A[r + c * n]. For each column j: copy row j (columns j+1..n-1, stride n)
// into column j (rows j+1..n-1, stride 1). Source and destination never overlap.
static void mirrorUpperToLower(cuSolverContext* ctx, double* d_A, int n) {
  for (int j = 0; j < n - 1; ++j) {
    CHECK_CUBLAS( cublasDcopy(ctx->cublasH, n - 1 - j,
                              d_A + j + static_cast<size_t>(j + 1) * n, n,
                              d_A + (j + 1) + static_cast<size_t>(j) * n, 1) );
  }
}

// QCidd = Q %*% Cidd on the GPU, where Cidd is the inverse currently stored (upper triangle) in ctx->d_C.
//
// Qt_in must be a dgCMatrix holding t(Q) (N x M). The CSC arrays of t(Q) are *exactly* the CSR arrays of Q
// (M x N): column pointers <-> row pointers, row indices <-> column indices. So Qt's p / i / x slots are uploaded
// as-is and described to cuSPARSE as a CSR matrix with M rows and N columns. No conversion or transposition needed.
//
// The number of stored entries of Q is taken from the first call and must stay the same afterwards (the values and
// the pattern may change; only the count is fixed).
//
// Returns a dgeMatrix (M x N, column-major), i.e. the same thing `Q %*% Cidd` gives on the CPU.
static Rcpp::S4 computeQCidd(cuSolverContext* ctx, SEXP Qt_in) {
  if (ctx->d_QCidd == nullptr) {
    Rcpp::stop("This context was created without Q: call cuSolverInit(N, rhs_cols_1, rhs_cols_2, M).");
  }
  if (!Rf_isS4(Qt_in)) {
    Rcpp::stop("Qt must be a dgCMatrix, i.e. t(Q) where Q is a dgCMatrix.");
  }
  Rcpp::S4 Qt(Qt_in);
  if (!Qt.is("dgCMatrix")) {
    Rcpp::stop("Qt must be a dgCMatrix, i.e. t(Q) where Q is a dgCMatrix.");
  }
  Rcpp::IntegerVector Qp   = Qt.slot("p");
  Rcpp::IntegerVector Qi   = Qt.slot("i");
  Rcpp::NumericVector Qx   = Qt.slot("x");
  Rcpp::IntegerVector Qdim = Qt.slot("Dim");
  
  const int N = ctx->N;     // dimension of Cidd = rows of Qt = columns of Q
  const int M = ctx->M;     // columns of Qt = rows of Q
  if (Qdim[0] != N || Qdim[1] != M) {
    Rcpp::stop("Qt is %i x %i but cuSolverInit was told to expect %i x %i.", Qdim[0], Qdim[1], N, M);
  }
  const int nnz = static_cast<int>(Qx.size());
  if (nnz == 0) {
    Rcpp::stop("Qt has no stored entries.");
  }
  
  if (ctx->nnzQ == 0) {
    // First call with Q: now the number of stored entries is known. Allocate the arrays that depend on it and create
    // the CSR descriptor of Q once; its pointers and nnz never change from here on.
    if (ctx->d_Qi) { cudaFree(ctx->d_Qi); ctx->d_Qi = nullptr; }
    if (ctx->d_Qx) { cudaFree(ctx->d_Qx); ctx->d_Qx = nullptr; }
    CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_Qi)), static_cast<size_t>(nnz) * sizeof(int)) );
    CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_Qx)), static_cast<size_t>(nnz) * sizeof(double)) );
    CHECK_CUSPARSE( cusparseCreateCsr(&(ctx->descA), M, N, nnz, ctx->d_Qp, ctx->d_Qi, ctx->d_Qx,
                                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) );
    ctx->nnzQ = nnz; // set last, so a failed first call leaves the context in the "not yet initialised for Q" state
  } else if (nnz != ctx->nnzQ) {
    Rcpp::stop("Qt has %i stored entries but the first call had %i; the number of non-zeros must stay the same.",
               nnz, ctx->nnzQ);
  }
  
  // Upload Q:
  CHECK_CUDA( cudaMemcpyAsync(ctx->d_Qp, Qp.begin(), (static_cast<size_t>(M) + 1) * sizeof(int), cudaMemcpyHostToDevice, ctx->stream) );
  CHECK_CUDA( cudaMemcpyAsync(ctx->d_Qi, Qi.begin(), static_cast<size_t>(nnz) * sizeof(int), cudaMemcpyHostToDevice, ctx->stream) );
  CHECK_CUDA( cudaMemcpyAsync(ctx->d_Qx, Qx.begin(), static_cast<size_t>(nnz) * sizeof(double), cudaMemcpyHostToDevice, ctx->stream) );
  
  // Fill the lower triangle of d_C so it holds the full symmetric Cidd.
  // (Phase 1 re-uploads C into d_C next iteration, and potri already destroyed the Cholesky factor, so nothing is lost.)
  mirrorUpperToLower(ctx, ctx->d_C, N);
  
  const double alpha = 1.0;
  const double beta  = 0.0;
  size_t bufSize = 0;
  CHECK_CUSPARSE( cusparseSpMM_bufferSize(ctx->cusparseH,
                                          CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          &alpha, ctx->descA, ctx->descB, &beta, ctx->descC,
                                          CUDA_R_64F, CUSPARSE_SPMM_ALG_DEFAULT, &bufSize) );
  ctx->d_spmmBuf.reserve(std::max<size_t>(bufSize, 16)); // grow-only, so this only allocates the first time (or if it needs more)
  CHECK_CUSPARSE( cusparseSpMM(ctx->cusparseH,
                               CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                               &alpha, ctx->descA, ctx->descB, &beta, ctx->descC,
                               CUDA_R_64F, CUSPARSE_SPMM_ALG_DEFAULT, ctx->d_spmmBuf.ptr) );
  
  // Copy the dense result back to R:
  const R_xlen_t len = static_cast<R_xlen_t>(M) * N;
  Rcpp::NumericVector QCidd_vec(Rcpp::no_init(len)); // no point zero-filling ~1 GB we are about to overwrite
  CHECK_CUDA( cudaMemcpyAsync(REAL(QCidd_vec), ctx->d_QCidd, len * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
  CHECK_CUDA( cudaStreamSynchronize(ctx->stream) );
  
  Rcpp::S4 QCidd_dge("dgeMatrix");
  QCidd_dge.slot("x") = QCidd_vec;
  QCidd_dge.slot("Dim") = Rcpp::IntegerVector::create(M, N);
  QCidd_dge.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  return QCidd_dge;
}

// [[Rcpp::export(cuSolverPhase1)]]
Rcpp::List cuSolverPhase1(SEXP ctx_ptr, SEXP C_in, SEXP RHS_in, bool returnCchol = false) {
  
  // Unwrap the persistent context:
  Rcpp::XPtr<cuSolverContext> ctx(ctx_ptr);
  if (ctx.get() == nullptr) Rcpp::stop("GPU pointer is null or has already been freed!");
  
  // Convert to armadillo matrices:
  arma::dmat C = convert_to_dmat(C_in);
  arma::dmat RHS = convert_to_dmat(RHS_in);
  
  // Left and ride side dimensions:
  int N = ctx->N;
  int rhs_cols = ctx->rhs_cols_1;
  if (static_cast<int>(RHS.n_cols) != rhs_cols) {
    Rcpp::stop("RHS dimension mismatch! Expected %i columns, got %i.", rhs_cols, RHS.n_cols);
  }
  
  // Assume upper storage ("U") from R's dsyMatrix/dtrMatrix conventions:
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
  
  // Push fresh data into the persistent buffers:
  CHECK_CUDA( cudaMemcpyAsync(ctx->d_C, C.memptr(), N * N * sizeof(double), cudaMemcpyHostToDevice, ctx->stream) );
  CHECK_CUDA( cudaMemcpyAsync(ctx->d_RHS, RHS.memptr(), N * rhs_cols * sizeof(double), cudaMemcpyHostToDevice, ctx->stream) );
  
  // In-place Cholesky decomposition (C becomes U or L):
  cusolverDnDpotrf(ctx->cusolverH, uplo, N, ctx->d_C, N, ctx->d_work, ctx->lwork, ctx->d_info);
  
  // Check if Cholesky factorization succeeded (info == 0):
  int info = 0;
  CHECK_CUDA( cudaMemcpyAsync(&info, ctx->d_info, sizeof(int), cudaMemcpyDeviceToHost, ctx->stream) );
  CHECK_CUDA( cudaStreamSynchronize(ctx->stream) ); // We must sync to ensure a potential error code is actually present in info!
  if (info != 0) {
    Rcpp::stop("The coefficient matrix C is not positive-definite (cuSOLVER info = %i).", info);
  }
  
  // Solve for bu in-place (cusolverDnDpotrs)
  // Overwrites d_RHS with the solution (bu):
  cusolverDnDpotrs(ctx->cusolverH, uplo, N, rhs_cols, ctx->d_C, N, ctx->d_RHS, N, ctx->d_info);
  
  // Copy bu back to the host:
  Rcpp::NumericVector bu_vec(N * rhs_cols);
  CHECK_CUDA( cudaMemcpyAsync(REAL(bu_vec), ctx->d_RHS, N * rhs_cols * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
  
  Rcpp::S4 bu_dge("dgeMatrix");
  bu_dge.slot("x") = bu_vec;
  bu_dge.slot("Dim") = Rcpp::IntegerVector::create(N, rhs_cols);
  bu_dge.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  
  // Compute log determinant: logDetC = 2 * sum(log(diag(C_factor)))
  cublasDcopy(ctx->cublasH, N, ctx->d_C, N + 1, ctx->d_diag, 1);
  
  // Copy only the N diagonal elements to the host:
  arma::dvec diag_host(N);
  CHECK_CUDA( cudaMemcpyAsync(diag_host.memptr(), ctx->d_diag, N * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
  CHECK_CUDA( cudaStreamSynchronize(ctx->stream) ); // We must sync to ensure diag_host actually holds the diagonal elements of the CF on the device!
  double logDetC = 0.0;
  for (int i = 0; i < N; ++i) {
    logDetC += std::log(diag_host[i]);
  }
  logDetC *= 2.0;
  
  // Copy the full Cholesky factor back to R only after convergence so we can make predictions and PEVs:
  // CHANGED: was a raw SEXP, which is not protected from R's garbage collector while List::create() allocates below.
  Rcpp::RObject chol_host = R_NilValue;
  if (returnCchol) {
    arma::dmat C_factor_host(N, N);
    CHECK_CUDA( cudaMemcpyAsync(C_factor_host.memptr(), ctx->d_C, N * N * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
    CHECK_CUDA( cudaStreamSynchronize(ctx->stream) ); // We must sync to ensure C_factor_host actually holds the elements from the CF on the device!
    C_factor_host = arma::trimatu(C_factor_host);
    chol_host = Rcpp::wrap(C_factor_host);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("bu") = bu_dge,
    // Rcpp::Named("z") = z_dge,
    Rcpp::Named("logDetC") = logDetC,
    Rcpp::Named("Cchol") = chol_host
  );
}

// CHANGED: optional third argument Qt_in. If it is a dgCMatrix holding t(Q), the returned list also contains
// QCidd = Q %*% Cidd (dgeMatrix, M x N); otherwise QCidd is NULL.
// [[Rcpp::export(cuSolverPhase2)]]
Rcpp::List cuSolverPhase2(SEXP ctx_ptr, SEXP XWjxZWj_in, SEXP Qt_in = R_NilValue) {
  
  // Unwrap the persistent context
  Rcpp::XPtr<cuSolverContext> ctx(ctx_ptr);
  if (ctx.get() == nullptr) Rcpp::stop("GPU pointer is null or has already been freed!");
  
  // Converting to an armadillo matrix:
  arma::dmat XWjxZWj = convert_to_dmat(XWjxZWj_in);
  
  // Left and right side dimensions:
  int N = ctx->N;
  int rhs_cols = ctx->rhs_cols_2;
  if (static_cast<int>(XWjxZWj.n_cols) != rhs_cols) {
    Rcpp::stop("XWjxZWj dimension mismatch! Expected %i columns, got %i.", rhs_cols, XWjxZWj.n_cols);
  }
  
  // Assume upper storage ("U") from R's dsyMatrix/dtrMatrix conventions:
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
  
  // Push the multiple right-hand side matrix (XWjxZWj) into the d_buWu buffer.
  // Note that XWjxZWj will be replaced by buWu:
  CHECK_CUDA( cudaMemcpyAsync(ctx->d_buWu, XWjxZWj.memptr(), N * rhs_cols * sizeof(double), cudaMemcpyHostToDevice, ctx->stream) );
  
  // Solve for buWu in-place using the existing Cholesky factor in d_C:
  cusolverDnDpotrs(ctx->cusolverH, uplo, N, rhs_cols, ctx->d_C, N, ctx->d_buWu, N, ctx->d_info);
  
  // Return buWu to R directly into a dgeMatrix
  Rcpp::NumericVector buWu_vec(N * rhs_cols);
  CHECK_CUDA( cudaMemcpyAsync(REAL(buWu_vec), ctx->d_buWu, N * rhs_cols * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
  
  Rcpp::S4 buWu_dge("dgeMatrix");
  buWu_dge.slot("x") = buWu_vec;
  buWu_dge.slot("Dim") = Rcpp::IntegerVector::create(N, rhs_cols);
  buWu_dge.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  
  // Invert C from its Cholesky factor:
  cusolverDnDpotri(ctx->cusolverH, uplo, N, ctx->d_C, N, ctx->d_work, ctx->lwork, ctx->d_info);
  
  // Allocate the raw R vector for the S4 dsyMatrix:
  Rcpp::NumericVector Ci_vec(N * N);
  
  // Beam C^-1 directly into the R vector
  CHECK_CUDA( cudaMemcpyAsync(REAL(Ci_vec), ctx->d_C, N * N * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
  
  // Sync here because the CPU cannot return the Ci_vec to R until the GPU finishes writing to it!
  CHECK_CUDA( cudaStreamSynchronize(ctx->stream) );
  
  // Construct the S4 dsyMatrix:
  Rcpp::S4 Ci_host("dpoMatrix");
  Ci_host.slot("x") = Ci_vec;
  Ci_host.slot("Dim") = Rcpp::IntegerVector::create(N, N);
  Ci_host.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  Ci_host.slot("uplo") = "U";
  Ci_host.slot("factors") = Rcpp::List::create(); // dpoMatrix strictly expects an empty list here by default
  
  // NEW: Q %*% Cidd via cuSPARSE. This must come after the Ci copy above, because it mirrors the upper triangle of d_C
  // into the lower triangle in place.
  Rcpp::RObject QCidd = R_NilValue;
  if (!Rf_isNull(Qt_in)) {
    QCidd = computeQCidd(ctx.get(), Qt_in);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("buWu") = buWu_dge,
    Rcpp::Named("Ci") = Ci_host,
    Rcpp::Named("QCidd") = QCidd
  );
}