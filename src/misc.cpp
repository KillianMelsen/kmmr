// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::S4 batchedWorkingVariates(const Rcpp::S4& umat_sigmai, const Rcpp::S4& batchedDerivs) {
  // Extract dimensions and data pointers:
  Rcpp::IntegerVector US_dim = umat_sigmai.slot("Dim");
  Rcpp::NumericVector US_x_ptr   = umat_sigmai.slot("x");
  Rcpp::IntegerVector D_dim = batchedDerivs.slot("Dim");
  Rcpp::NumericVector D_x_ptr   = batchedDerivs.slot("x");
  
  // Wrap data into arma::mat without copying:
  arma::mat US(US_x_ptr.begin(), US_dim[0], US_dim[1], false, true);
  arma::mat D(D_x_ptr.begin(), D_dim[0], D_dim[1], false, true);
  
  // Dimensions:
  const arma::uword g = US.n_rows; // Dimension of non-variable structue, e.g., genotypes.
  const arma::uword p = US.n_cols; // Dimension of variable structure, e.g., environments.
  const arma::uword q = D.n_cols / p; // Number of derivatives batched together in D.
  
  // Batched execution:
  Rcpp::NumericVector out_x(g * p * q); // Pre-allocate memory for the result.
  arma::mat Y(out_x.begin(), g, p * q, false, true); // Wrap in arma::mat without copying.
  Y = US * D;
  
  // Return as dgeMatrix for multiplying with design matrix:
  Rcpp::S4 out("dgeMatrix");
  out.slot("x") = out_x;
  out.slot("Dim") = Rcpp::IntegerVector::create(g * p, q);
  out.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector batchedScores(const Rcpp::S4& basisScores, const Rcpp::S4& batchedDerivs) {
  // Extract dimensions and data pointers:
  Rcpp::IntegerVector D_dim = batchedDerivs.slot("Dim");
  const arma::uword p = D_dim[0];
  const arma::uword q  = D_dim[1] / p;
  Rcpp::NumericVector D_x_ptr = batchedDerivs.slot("x");
  Rcpp::NumericVector B_x_ptr = basisScores.slot("x");
  
  // Wrap data into arma::mat without copying, and reshape simultaneously to flatten:
  arma::mat D_flat(D_x_ptr.begin(), p * p, q, false, true); // each column represents an entire derivative.
  arma::mat B(B_x_ptr.begin(), p, p, false, true); // First 
  
  // basisScores is a dsyMatrix which contains nonsense in one of the triangles, so we symmetrize in-place:
  std::string uplo = Rcpp::as<std::string>(basisScores.slot("uplo"));
  if (uplo == "U") { // Upper triangle is valid. Copy upper to lower.
    for (arma::uword c = 0; c < p; ++c) {
      for (arma::uword r = c + 1; r < p; ++r) {
        B(r, c) = B(c, r);
      }
    }
  } else { // Lower triangle is valid. Copy lower to upper.
    for (arma::uword c = 0; c < p; ++c) {
      for (arma::uword r = c + 1; r < p; ++r) {
        B(c, r) = B(r, c);
      }
    }
  }
  
  // Now flatten into a vector that contains all basis scores:
  arma::vec b_flat(B_x_ptr.begin(), p * p, false, true);
  
  // Compute all scores simultaneously as a single matrix - vector multiplication.
  // Each score is the inner product of a vector representing n entire derivative and a vector representing all basis scores.
  // Returning as a base R vector is fine as we'll put them into the scores dgeMatrix manually:
  return Rcpp::wrap(D_flat.t() * b_flat);
}


