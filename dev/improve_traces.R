# Summing takes about 0.125 seconds total
# Sub-setting takes about 3.2 seconds total
# We need nX, sigmai, dim.nvar, nvarGi, C
rm(list = ls())
load("dev/improve_traces.RData")
re <- 1
corner <- nX + 1

# Option 1, 3.2 seconds ====
# traces <- Matrix(0, nrow(sigmai[[re]]), ncol(sigmai[[re]]), sparse = FALSE, doDiag = FALSE)
# system.time({
# for (i in 1:nrow(sigmai[[re]])) {
#   for (j in 1:ncol(sigmai[[re]])) {
#     if (j >= i) { # Only do the upper triangle
#       rstart <- (corner + (i - 1) * dim.nvar[re])
#       rend <- rstart - 1 + dim.nvar[re]
#       cstart <- (corner + (j - 1) * dim.nvar[re])
#       cend <- cstart - 1 + dim.nvar[re]
#       # I suspect this can be done faster for a sparse inverse subset (C) by selecting the relevant (nonzero) elements from the current block of C
#       # and the corresonding elements of nvarGi. But for now let's keep the simple subsetting and hadamard product.
#       traces[i, j] <- sum(nvarGi[[re]] * C[rstart:rend, cstart:cend])
#     }
#   }
# }
# })

library(Matrix)
library(tictoc)

# Option 2, 0.5 seconds ====
# # 1. Define dimensions and extract the relevant sub-matrix ONE time
# system.time({
# p <- nrow(sigmai[[re]])
# idx_start <- corner
# idx_end <- corner + p * dim.nvar - 1
# 
# # Subsetting C once is vastly faster than doing it p^2 times
# tic("Subsetting")
# C_sub <- C[idx_start:idx_end, idx_start:idx_end, drop = FALSE]
# toc()
# 
# # 2. Coerce to a general triplet matrix (dgTMatrix)
# # We coerce to a general matrix first so BOTH triangles are explicitly represented.
# # This prevents missing off-diagonal non-zeros that might be stored in the lower triangle.
# tic("Converting C_sub")
# C_T <- as(as(C_sub, "generalMatrix"), "TsparseMatrix")
# toc()
# 
# # 3. Calculate the block indices and within-block indices
# # IMPORTANT: @i and @j in TsparseMatrix are 0-based!
# tic("Computing indices")
# block_i <- (C_T@i %/% dim.nvar) + 1
# block_j <- (C_T@j %/% dim.nvar) + 1
# 
# in_block_i <- (C_T@i %% dim.nvar) + 1
# in_block_j <- (C_T@j %% dim.nvar) + 1
# toc()
# 
# # 4. Extract corresponding elements from nvarGi and multiply
# # We subset the dense matrix nvarGi using a 2-column matrix of coordinates
# tic("Extracting nvarGi values")
# G_vals <- nvarGi[[re]][cbind(in_block_i, in_block_j)]
# toc()
# tic("Computing values")
# vals <- C_T@x * G_vals
# toc()
# 
# # 5. Filter for the upper triangle (j >= i) to match your original logic
# tic("Filtering")
# valid <- block_j >= block_i
# toc()
# 
# # 6. Sum the values into the traces matrix using sparseMatrix aggregation
# # A beautiful trick in R: `sparseMatrix` automatically sums values 
# # if there are duplicate (i, j) coordinate pairs!
# tic("Creating traces matrix")
# traces_sparse <- sparseMatrix(
#   i = block_i[valid],
#   j = block_j[valid],
#   x = vals[valid],
#   dims = c(p, p)
# )
# toc()
# 
# # Convert back to a standard dense matrix
# tic("Conversion of traces matrix to dgeMatrix")
# traces <- as(traces_sparse, "dgeMatrix")
# toc()
# })
# all.equal(traces, correctTraces)

# Option 3 0.1-0.15 seconds ====
# system.time({
# p <- nrow(sigmai[[re]])
# d <- as.integer(nrow(nvarGi[[re]])) 
# idx_start <- corner
# idx_end <- corner + p * d - 1
# 
# # C_sub is a dscMatrix (lower triangular symmetric)
# C_sub <- C[idx_start:idx_end, idx_start:idx_end, drop = FALSE]
# 
# # Reconstruct coordinates
# i_vec <- C_sub@i
# j_vec <- rep(0L:(ncol(C_sub) - 1L), diff(C_sub@p))
# 
# block_i <- (i_vec %/% d) + 1L
# block_j <- (j_vec %/% d) + 1L
# 
# in_block_i <- i_vec %% d
# in_block_j <- j_vec %% d
# 
# # Convert nvarGi to a standard matrix for safe, fast 1D indexing
# G_mat <- as.matrix(nvarGi[[re]])
# idx_1d <- in_block_i + (in_block_j * d) + 1L
# 
# # Symmetry trick: double the off-diagonal elements inside the diagonal blocks
# mult <- rep(1.0, length(i_vec))
# mult[block_i == block_j & in_block_i > in_block_j] <- 2.0
# 
# vals <- C_sub@x * G_mat[idx_1d] * mult
# 
# # Note: We assign row = block_j, col = block_i. 
# # Because C_sub is lower triangular (block_i >= block_j), 
# # this perfectly populates the UPPER triangle of traces!
# traces_sparse <- sparseMatrix(
#   i = block_j,
#   j = block_i,
#   x = vals,
#   dims = c(p, p)
# )
# traces <- as(traces_sparse, "dgeMatrix")
# })
# all.equal(correctTraces, traces)

# Option 4, x seconds ====
library(Rcpp)
cppFunction('
NumericMatrix get_traces_sym(IntegerVector p_ptr, IntegerVector i_ptr, NumericVector x_ptr,
                             NumericMatrix G, int p, int d) {

    NumericMatrix traces(p, p);

    for (int col = 0; col < p_ptr.length() - 1; ++col) {
        int block_j = col / d;
        int in_j = col % d;

        int start = p_ptr[col];
        int end = p_ptr[col + 1];

        for (int idx = start; idx < end; ++idx) {
            int row = i_ptr[idx];
            int block_i = row / d;
            int in_i = row % d;

            double val = x_ptr[idx] * G(in_i, in_j);

            // If we are calculating a block on the diagonal...
            if (block_i == block_j) {
                // ... double the lower off-diagonal elements to account for the missing upper half
                if (in_i > in_j) {
                    val *= 2.0;
                }
            }

            // Assign to (block_j, block_i) to form the upper triangle
            traces(block_j, block_i) += val;
        }
    }
    return traces;
}
')

# --- Usage ---
p <- nrow(sigmai[[re]])
d <- nrow(nvarGi[[re]])
idx_start <- corner
idx_end <- corner + p * d - 1

system.time({
C_sub <- C[idx_start:idx_end, idx_start:idx_end, drop = FALSE]

# Note: wrap nvarGi in as.matrix() to ensure it strips the `dsymatrix` class
# so Rcpp can interpret it natively as a NumericMatrix.
traces <- get_traces_sym(
  p_ptr = C_sub@p,
  i_ptr = C_sub@i,
  x_ptr = C_sub@x,
  G = as.matrix(nvarGi[[re]]),
  p = p,
  d = d
)
})
all.equal(forceSymmetric(traces), correctTraces)

