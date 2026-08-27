// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

// Device properties function
//[[Rcpp::export(.getDeviceProps)]]
Rcpp::List getDeviceProps(const int& id) {
  
  // Create CUDA device properties variable
  cudaDeviceProp prop;
  
  // Query GPU properties
  cudaGetDeviceProperties(&prop, id);
  
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
      
      // 1. Initialize dmat instantly using the memory pointer
      // (This safely copies the data by default, preventing accidental modification of the R object)
      arma::dmat out(data.begin(), dims[0], dims[1]);
      
      // 2. Symmetrize if it is a dsyMatrix
      if (obj.is("dsyMatrix") && obj.hasSlot("uplo")) {
        std::string uplo = Rcpp::as<std::string>(obj.slot("uplo"));
        if (uplo == "U") {
          out = arma::symmatu(out); // Copies Upper to Lower
        } else if (uplo == "L") {
          out = arma::symmatl(out); // Copies Lower to Upper
        }
      }
      
      return out;
    }
  }
  Rcpp::stop("Input must be a dense Matrix object (dgeMatrix or dsyMatrix)...");
}

// Helper function to extract and cast to fmat
arma::fmat convert_to_fmat(SEXP x) {
  if (Rf_isS4(x)) {
    Rcpp::S4 obj(x);
    if (obj.hasSlot("x") && obj.hasSlot("Dim")) {
      Rcpp::NumericVector data = obj.slot("x");
      Rcpp::IntegerVector dims = obj.slot("Dim");
      
      // 1. Initialize fmat and cast double to float in one pass
      arma::fmat out(dims[0], dims[1]);
      for(int i = 0; i < data.size(); ++i) {
        out[i] = static_cast<float>(data[i]);
      }
      
      // 2. Symmetrize if it is a dsyMatrix
      if (obj.is("dsyMatrix") && obj.hasSlot("uplo")) {
        std::string uplo = Rcpp::as<std::string>(obj.slot("uplo"));
        if (uplo == "U") {
          out = arma::symmatu(out); // Copies Upper to Lower
        } else if (uplo == "L") {
          out = arma::symmatl(out); // Copies Lower to Upper
        }
      }
      
      return out;
    }
  }
  Rcpp::stop("Input must be a dense Matrix object (dgeMatrix or dsyMatrix)...");
}

void gc_free_gpu_memory(double* ptr) {
  if (ptr != nullptr) cudaFree(ptr);
}

// [[Rcpp::export(gpuFreePointer)]]
void gpuFreePointer(SEXP ptr_in) {
  Rcpp::XPtr<double> ptr(ptr_in);
  double* d_C = ptr.get();
  if (d_C != nullptr) {
    cudaFree(d_C);
    R_ClearExternalPtr(ptr_in); // Crucial: Stops R's GC from double-freeing later
  }
}

// [[Rcpp::export(cuSolverPhase1)]]
Rcpp::List cuSolverPhase1(SEXP C_in, SEXP RHS_in, bool returnCchol = false) {
  
  // Convert to armadillo matrices:
  arma::dmat C = convert_to_dmat(C_in);
  arma::dmat RHS = convert_to_dmat(RHS_in);
  
  // Left and ride side dimensions:
  int N = C.n_rows;
  int rhs_cols = RHS.n_cols;
  
  // Create handles and error variables:
  cudaError_t err;
  cublasHandle_t cublasH = NULL;
  cusolverDnHandle_t cusolverH = NULL;
  cudaStream_t stream = NULL;
  
  cublasCreate(&cublasH);
  cusolverDnCreate(&cusolverH);
  cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  cublasSetStream(cublasH, stream);
  cusolverDnSetStream(cusolverH, stream);
  
  // Allocate device memory for C and RHS, then copy from host:
  double *d_C = nullptr;
  double *d_RHS = nullptr;
  int *d_info = nullptr; // To check if matrix is positive-definite
  
  cudaMalloc(reinterpret_cast<void**>(&d_C), N * N * sizeof(double));
  cudaMalloc(reinterpret_cast<void**>(&d_RHS), N * rhs_cols * sizeof(double));
  cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int));
  
  cudaMemcpyAsync(d_C, C.memptr(), N * N * sizeof(double), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(d_RHS, RHS.memptr(), N * rhs_cols * sizeof(double), cudaMemcpyHostToDevice, stream);
  
  // Cholesky Factorization (cusolverDnDpotrf)
  // Assume upper storage ("U") from R's dsyMatrix/dtrMatrix conventions:
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
  
  int lwork = 0;
  cusolverDnDpotrf_bufferSize(cusolverH, uplo, N, d_C, N, &lwork);
  
  double *d_work = nullptr;
  cudaMalloc(reinterpret_cast<void**>(&d_work), lwork * sizeof(double));
  
  // In-place Cholesky decomposition (C becomes U or L):
  cusolverDnDpotrf(cusolverH, uplo, N, d_C, N, d_work, lwork, d_info);
  
  // Check if Cholesky factorization succeeded (info == 0):
  int info = 0;
  cudaMemcpyAsync(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost, stream);
  cudaStreamSynchronize(stream); // We must sync to ensure a potential error code is actually present in info!
  if (info != 0) {
    cudaFree(d_C);
    cudaFree(d_RHS);
    cudaFree(d_info);
    cudaFree(d_work);
    cublasDestroy(cublasH);
    cusolverDnDestroy(cusolverH);
    cudaStreamDestroy(stream);
    Rcpp::stop("The coefficient matrix C is not positive-definite (cuSOLVER info = %i).", info);
  }
  cudaFree(d_work);
  
  // Solve for bu in-place (cusolverDnDpotrs)
  // Overwrites d_RHS with the solution (bu):
  cusolverDnDpotrs(cusolverH, uplo, N, rhs_cols, d_C, N, d_RHS, N, d_info);
  
  // Copy bu back to the host:
  arma::dmat bu_host(N, rhs_cols);
  cudaMemcpyAsync(bu_host.memptr(), d_RHS, N * rhs_cols * sizeof(double), cudaMemcpyDeviceToHost, stream);
  
  // Compute z = U * bu in-place using cuBLAS triangular matrix vector multiply (cublasDtrmv)
  // Note: z needs to be computed using the Cholesky factor still sitting in d_C and the solved bu in d_RHS.
  cublasDtrmv(cublasH, uplo, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, N, d_C, N, d_RHS, 1);
  
  // Copy z back to the host:
  arma::dmat z_host(N, rhs_cols);
  cudaMemcpyAsync(z_host.memptr(), d_RHS, N * rhs_cols * sizeof(double), cudaMemcpyDeviceToHost, stream);
  
  // Compute log determinant: logDetC = 2 * sum(log(diag(C_factor)))
  double *d_diag = nullptr;
  cudaMalloc(reinterpret_cast<void**>(&d_diag), N * sizeof(double));
  cublasDcopy(cublasH, N, d_C, N + 1, d_diag, 1);
  
  // Copy only the N diagonal elements to the host:
  arma::dvec diag_host(N);
  cudaMemcpyAsync(diag_host.memptr(), d_diag, N * sizeof(double), cudaMemcpyDeviceToHost, stream);
  cudaStreamSynchronize(stream); // We must sync to ensure diag_host actually holds the diagonal elements of the CF on the device!
  double logDetC = 0.0;
  for (int i = 0; i < N; ++i) {
    logDetC += std::log(diag_host[i]);
  }
  logDetC *= 2.0;
  cudaFree(d_diag);
  
  // Copy the full Cholesky factor back to R only after convergence so we can make predictions and PEVs:
  SEXP chol_host = R_NilValue;
  if (returnCchol) {
    arma::dmat C_factor_host(N, N);
    cudaMemcpyAsync(C_factor_host.memptr(), d_C, N * N * sizeof(double), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream); // We must sync to ensure C_factor_host actually holds the elements from the CF on the device!
    if (uplo == CUBLAS_FILL_MODE_UPPER) {
      C_factor_host = arma::trimatu(C_factor_host);
    } else {
      C_factor_host = arma::trimatl(C_factor_host);
    }
    chol_host = Rcpp::wrap(C_factor_host);
  }
  
  // Cleanup (keep d_C alive inside XPtr!):
  cudaFree(d_RHS);
  cudaFree(d_info);
  cublasDestroy(cublasH);
  cusolverDnDestroy(cusolverH);
  cudaStreamDestroy(stream);
  
  // Wrap d_C in XPtr for phase 2:
  Rcpp::XPtr<double> ptr_out(d_C, true, gc_free_gpu_memory);
  
  return Rcpp::List::create(
    Rcpp::Named("bu") = bu_host,
    Rcpp::Named("z") = z_host,
    Rcpp::Named("logDetC") = logDetC,
    Rcpp::Named("Cchol") = chol_host,
    Rcpp::Named("ptr") = ptr_out
  );
}










