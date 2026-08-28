// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#define CHECK_CUDA(call) {                                                                          \
  cudaError_t err = call;                                                                           \
  if (err != cudaSuccess) {                                                                         \
    Rcpp::stop("CUDA error at line %d of /src/gpuCode.cpp: %s", __LINE__, cudaGetErrorString(err)); \
  }                                                                                                 \
}                                                                                                 

// The Persistent Context Structure:
struct cuSolverContext {
  cublasHandle_t cublasH;
  cusolverDnHandle_t cusolverH;
  cudaStream_t stream;
  
  int N;
  int rhs_cols;
  int lwork; // Stores the max workspace size needed
  
  // Persistent memory blocks
  double *d_C;
  double *d_RHS;
  double *d_buWu;
  double *d_work;
  double *d_diag;
  int *d_info;
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
      if (ctx->cublasH) cublasDestroy(ctx->cublasH);
      if (ctx->cusolverH) cusolverDnDestroy(ctx->cusolverH);
      if (ctx->stream) cudaStreamDestroy(ctx->stream);
      
      delete ctx; // Free the struct itself
      R_ClearExternalPtr(ptr_in); // Prevent double-free
    }
  }
}

// The initializer (called from R once before iterating):
// [[Rcpp::export(cuSolverInit)]]
SEXP cuSolverInit(int N, int rhs_cols) {
  cuSolverContext* ctx = new cuSolverContext();
  ctx->N = N;
  ctx->rhs_cols = rhs_cols;
  
  cublasCreate(&(ctx->cublasH));
  cusolverDnCreate(&(ctx->cusolverH));
  CHECK_CUDA( cudaStreamCreateWithFlags(&(ctx->stream), cudaStreamNonBlocking) );
  cublasSetStream(ctx->cublasH, ctx->stream);
  cusolverDnSetStream(ctx->cusolverH, ctx->stream);
  
  // Allocate all memory upfront:
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_C)), N * N * sizeof(double)) );
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_RHS)), N * rhs_cols * sizeof(double)) );
  CHECK_CUDA( cudaMalloc(reinterpret_cast<void**>(&(ctx->d_buWu)), N * rhs_cols * sizeof(double)) );
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

// Helper function to extract and cast to fmat:
arma::fmat convert_to_fmat(SEXP x) {
  if (Rf_isS4(x)) {
    Rcpp::S4 obj(x);
    if (obj.hasSlot("x") && obj.hasSlot("Dim")) {
      Rcpp::NumericVector data = obj.slot("x");
      Rcpp::IntegerVector dims = obj.slot("Dim");
      
      // 1. Allocate the float matrix
      arma::fmat out(dims[0], dims[1]);
      
      // 2. Grab raw pointers for maximum compiler vectorization speed
      const double* src = data.begin();
      float* dest = out.memptr();
      arma::uword total_elements = out.n_elem;
      
      // 3. Blazing fast cast loop
      for(arma::uword i = 0; i < total_elements; ++i) {
        dest[i] = static_cast<float>(src[i]);
      }
      
      return out;
    }
  }
  Rcpp::stop("Input must be a dense Matrix object (dgeMatrix or dsyMatrix)...");
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
  int rhs_cols = ctx->rhs_cols;
  
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
  
  // Compute z = U * bu in-place using cuBLAS triangular matrix vector multiply (cublasDtrmv)
  // Note: z needs to be computed using the Cholesky factor still sitting in d_C and the solved bu in d_RHS.
  cublasDtrmv(ctx->cublasH, uplo, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, N, ctx->d_C, N, ctx->d_RHS, 1);
  
  // Copy z back to the host:
  Rcpp::NumericVector z_vec(N * rhs_cols);
  CHECK_CUDA( cudaMemcpyAsync(REAL(z_vec), ctx->d_RHS, N * rhs_cols * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
  
  Rcpp::S4 z_dge("dgeMatrix");
  z_dge.slot("x") = z_vec;
  z_dge.slot("Dim") = Rcpp::IntegerVector::create(N, rhs_cols);
  z_dge.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  
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
  SEXP chol_host = R_NilValue;
  if (returnCchol) {
    arma::dmat C_factor_host(N, N);
    CHECK_CUDA( cudaMemcpyAsync(C_factor_host.memptr(), ctx->d_C, N * N * sizeof(double), cudaMemcpyDeviceToHost, ctx->stream) );
    CHECK_CUDA( cudaStreamSynchronize(ctx->stream) ); // We must sync to ensure C_factor_host actually holds the elements from the CF on the device!
    C_factor_host = arma::trimatu(C_factor_host);
    chol_host = Rcpp::wrap(C_factor_host);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("bu") = bu_dge,
    Rcpp::Named("z") = z_dge,
    Rcpp::Named("logDetC") = logDetC,
    Rcpp::Named("Cchol") = chol_host
  );
}










