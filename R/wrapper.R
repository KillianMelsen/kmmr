#' RotateInitL
#'
#' @param L Matrix of loadings
#'
#' @returns The rotated loadings matrix so that the upper triangle is zero.
#' @keywords internal
rotateInitL <- function(L) {
  L.rotated <- L
  m <- ncol(L.rotated)
  p <- nrow(L.rotated)
  coords <- which(upper.tri(matrix(0, p, m)), arr.ind = TRUE)
  givens <- vector("list", nrow(coords))
  for (i in 1:nrow(coords)) {
    row <- coords[i, "row"]
    col <- coords[i, "col"]
    givens[[i]] <- diag(m)
    h <- sqrt(sum(c(L.rotated[row, row], L.rotated[row, col])^2))
    c <- L.rotated[row, row] / h
    s <- L.rotated[row, col] / h
    diag(givens[[i]])[c(row, col)] <- c
    givens[[i]][row, col] <- -s
    givens[[i]][col, row] <- s
    L.rotated <- L.rotated %*% givens[[i]]
  }
  L.rotated[upper.tri(L.rotated)] <- 0 # Make sure the upper triangle is exactly 0 and not xe-18
  rotmat <- Reduce("%*%", givens)
  return(L.rotated)
}

#' revInteractions
#'
#' @param form Formula.
#'
#' @returns The formula with kronecker products/interaction terms reversed.
#' @keywords internal
revInteractions <- function(form) {
  tt <- terms(form)
  response <- deparse(form[[2]])
  term_labels <- attr(tt, "term.labels")
  has_intercept <- attr(tt, "intercept") == 1
  
  term_labels <- vapply(term_labels, function(term) {
    if (grepl(":", term, fixed = TRUE)) {
      paste(rev(strsplit(term, ":", fixed = TRUE)[[1]]), collapse = ":")
    } else {
      term
    }
  }, character(1))
  
  rhs <- paste(c(ifelse(has_intercept, "1", "-1"), term_labels),
               collapse = " + ")
  
  return(as.formula(paste(response, "~", rhs)))
}

#' re
#'
#' @param formula Random effect formula.
#' @param vars To which factors should the variances be applied?
#' @param ... 
#'
#' @returns A list with all information on the structure of the random effect, to be used with the solver.
#' @keywords internal
re <- function(formula, vars = ~ 1, ...) {
  data <- list(...)$data
  calls <- strsplit(as.character(formula[2]), ":")[[1]] # vector of function calls for each structure/factor.
  results <- lapply(calls, function(call) {
    callStack <- grepl("stack\\(.*", call)
    if (callStack) {
      expr <- parse(text = call)[[1]]
      expr$data <- data
      eval(expr)
    } else {
      eval(parse(text = call), envir = data)
    }
  })
  isStack <- unlist(lapply(results, function(element) {"sizeOfStack" %in% names(element)}))
  if (any(isStack)) {
    if (sum(isStack) > 1) {
      stop("Only one stack per random effect allowed.")
    }
    if (!(all(isStack == c(TRUE, rep(FALSE, length(isStack) - 1))))) {
      stop("The stack must come first in the random effect specification.")
    }
    sos <- results[[which(isStack)]]$sizeOfStack
  }
  
  if (any(isStack)) {
    struc <- c(list(gsub("(.*)\\(.*$", "\\1", strsplit(gsub("stack\\(\\s*~\\s*(.*)\\)", "\\1", calls[isStack]), "\\s*\\+\\s*")[[1]])),
               gsub("(.*)\\(.*$", "\\1", calls[!isStack]))
  } else {
    struc <- gsub("(.*)\\(.*$", "\\1", calls) # What are the correlation structures for this factor?
  }
  factorNames <- gsub(".*\\(([^,\\)]+).*$", "\\1", calls) # Which factors are we evaluating?
  
  bases <- lapply(results, `[[`, "base")
  levels <- lapply(results, `[[`, "levels")
  nlevels <- vapply(results, `[[`, numeric(1), "nlevels")
  # nparams <- vapply(results, `[[`, numeric(1), "nparams")
  nparams <- lapply(results, `[[`, "nparams")
  part <- vapply(results, `[[`, character(1), "part")
  # types <- unlist(lapply(results, `[[`, "types"))
  types <- lapply(results, `[[`, "types")
  # indexCor <- unlist(lapply(results, `[[`, "indexCor"), recursive = FALSE) # I don't think we ever needed this recursive = FALSE...
  indexCor <- lapply(results, `[[`, "indexCor")
  forFactor <- lapply(seq_along(results), function(i) {
    if (isStack[i]) {
      lapply(results[[i]]$nparams, function(x) rep(i, x))
    } else {
      rep(i, results[[i]]$nparams)
    }
    })
  if (any(isStack)) {
    # Parameters that are part of a stack have an associated stackLayer number
    stackLayer <- c(unlist(lapply(seq_along(nparams[[which(isStack)]]), function(i) {rep(i, nparams[[which(isStack)]][[i]])})),
                    unlist(lapply(nparams[which(!isStack)], function(x) {rep(0, x)})))
  } else {
    # 0 Means the parameter is not part of a stack
    stackLayer <- unlist(lapply(nparams, function(x) {rep(0, x)}))
  }
  
  if (as.character(vars[2]) == "1") {
    applyvars <- 0
    if (any(isStack)) {
      nvars <- sos
    } else {
      nvars <- 1
    }
  } else {
    applyvars <- which(factorNames %in% strsplit(as.character(vars[2]), "[:]")[[1]])
    if (any(isStack)) {
      nvars <- prod(nlevels[applyvars]) * sos
    } else {
      nvars <- prod(nlevels[applyvars])
    }
  }
  
  # Determine which part of the random effect is variable and which part is non-variable:
  i.var <- setdiff(unique(c(applyvars, which(part == "var"))), 0)
  i.nvar <- setdiff(1:length(calls), i.var)
  if (!(all(i.var == 1:length(i.var)))) {
    stop(paste0("In '~ re(~ ", as.character(formula[2]), ", vars = ~ ", as.character(vars[2]), ")', linear kernels to which no variances are applied must be placed in the right part of the kronecker product."))
  }
  
  # Now we can add the variances to all vectors:
  types <- c(rep("V", nvars), na.omit(unlist(types)))
  if (any(isStack)) {
    stackLayer <- c(rep(1:sos, each = nvars / sos), stackLayer)
  } else {
    stackLayer <- c(rep(0, nvars), stackLayer)
  }
  forFactor <- c(rep(0, nvars), unlist(forFactor))
  if (any(isStack)) {
    nparams[[which(isStack)]] <- unlist(nparams[[which(isStack)]])
  }
  nparams.list <- c(list(nvars), nparams)
  # nparams <- c(nvars, unlist(lapply(seq_along(nparams), function(i) {sum(nparams[[i]])})))
  nparams <- unlist(lapply(nparams.list, sum))
  
  # Possible option for indexCor handling:
  # We need to unlist exactly twice to handle corg() and cors() structures appropriately.
  # There is probably a more elegant way to do this...
  # This is because corg() returns indexCor as a vector as each parameter is associated
  # with exactly one correlation.
  # cors() on the other hand returns indexCor as a list containing a  vector because the single
  # parameter is associated with all correlations.
  # If we unlist exactly twice we will get a single vector for corg(), and a list of vectors
  # for cors.
  # This ended up not working for something like re(~ cors(Env):lk(Gen, Ki), vars = ~ 1), I suspect
  # because if there is no stack, unlisting twice means we are actually unlisting the cors indexCor list into a vector like for corg.
  # I think a simple if statement will solve it.
  if (any(isStack)) { # If there is a stack we need to unlist exactly twice.
    indexCor <- unlist(unlist(indexCor, recursive = FALSE), recursive = FALSE)
  } else { # If there is no stack we need to unlist only once to avoid unlisting the cors list(1:p)
    indexCor <- unlist(indexCor, recursive = FALSE)
  }
  # Here we omit NAs (for liner kernels with no associated parameters for example)
  indexCor <- indexCor[!is.na(indexCor)] # [!is.na(indexCor)] works for both a list (cors) and a vector (corg)
  indexCor <- c(rep(0, nvars), indexCor) # And finally adding the 0s for the variances. This also works for both lists and vectors.
  
  indices <- vector("list", 1 + length(calls))
  index <- 1
  for (i in 1:(1 + length(calls))) {
    if (nparams[i] > 0) {
      indices[[i]] <- index:(index - 1 + nparams[i])
      index <- index + nparams[i]
    }
  }
  
  # Construct the design matrix:
  superlevels <- do.call(paste, c(rev(expand.grid(rev(levels))), sep = ":"))
  superfactor <- factor(do.call(paste, c(data[factorNames], sep = ":")), superlevels)
  Z <- sparse.model.matrix(~ -1 + superfactor)
  colnames(Z) <- superlevels
  params <- character(length(types))
  
  # Figure out the labels for the variances:
  h <- 1
  if (applyvars[1] == 0) {
    if (any(isStack)) {
      params[1:nvars] <- sprintf("vars{1|stack-%d}", 1:sos)
      h <- h + nvars
    } else {
      params[1] <- "vars{1}"
      h <- h + 1
    }
  } else {
    labels <- do.call(paste, c(rev(expand.grid(rev(levels[applyvars]))), sep = ":"))
    p <- length(labels)
    if (any(isStack)) {
      for (s in 1:sos) {
        params[h:(h - 1 + p)] <- sprintf("vars{%s-%s|stack-%d}", paste(factorNames[applyvars], collapse = ":"), labels, s)
        h <- h + p
      }
    } else {
      params[1:p] <- sprintf("vars{%s-%s}", paste(factorNames[applyvars], collapse = ":"), labels)
      h <- h + p
    }
  }
  
  # Figure out the labels for the other factors/structures:
  for (i in 1:length(factorNames)) {
    if (nparams[i + 1] > 0) { # + 1 because we can skip the variances
      if (isStack[i]) {
        for (s in 1:sos) {
          p <- nparams.list[[i + 1]][s]
          structure <- struc[[i]][s]
          if (structure == "corg") {
            labels <- t(outer(levels[[i]], levels[[i]], function(x, y) paste(x, y, sep = ":")))[lower.tri(matrix(1, nlevels[i], nlevels[i]))]
            params[h:(h - 1 + p)] <- sprintf("corg{%s-%s|stack-%d}", factorNames[i], labels, s)
            h <- h + p
          } else if (structure == "cors") {
            params[h] <- sprintf("cors{%s|stack-%d}", factorNames[i], s)
            h <- h + 1
          } else if (structure == "gk") {
            params[h] <- sprintf("bw{%s|stack-%d}", factorNames[i], s)
            h <- h + 1
          } else if (structure %in% c("fa", "exfa")) {
            m <- p / nlevels[i]
            labels <- paste(rep(levels[[i]], m), rep(sprintf("F%d", 1:m), each = nlevels[i]), sep = ":")
            params[h:(h - 1 + p)] <- sprintf("loading{%s-%s|stack-%d}", factorNames[i], labels, s)
            h <- h + p
          }
        }
      } else {
        p <- nparams[i + 1]
        structure <- struc[i]
        if (structure == "corg") {
          labels <- t(outer(levels[[i]], levels[[i]], function(x, y) paste(x, y, sep = ":")))[lower.tri(matrix(1, nlevels[i], nlevels[i]))]
          params[h:(h - 1 + p)] <- sprintf("corg{%s-%s}", factorNames[i], labels)
          h <- h + p
        } else if (structure == "cors") {
          params[h] <- sprintf("cors{%s}", factorNames[i])
          h <- h + 1
        } else if (structure == "gk") {
          params[h] <- sprintf("bw{%s}", factorNames[i])
          h <- h + 1
        } else if (structure %in% c("fa", "exfa")) {
          m <- p / nlevels[i]
          labels <- paste(rep(levels[[i]], m), rep(sprintf("F%d", 1:m), each = nlevels[i]), sep = ":")
          params[h:(h - 1 + p)] <- sprintf("loading{%s-%s}", factorNames[i], labels)
          h <- h + p
        }
      }
    }
  }
  
  return(list(Z = Z,
              indices = indices,
              applyvars = applyvars,
              types = types,
              nlevels = nlevels,
              forFactor = forFactor,
              indexCor = indexCor,
              bases = bases,
              struc = struc,
              isStack = isStack,
              nparams = sum(nparams),
              nparams.list = nparams.list,
              i.var = i.var,
              i.nvar = i.nvar,
              params = params))
}

#' \code{corg} is for a general correlation structure
#'
#' @param factor Data factor.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
corg <- function(factor) {
  levels <- levels(factor)
  p <- length(levels)
  nparams <- (p^2 - p) / 2
  return(list(base = NA,
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = rep("R", nparams),
              indexCor = 1:nparams,
              part = "var"))
}

#' \code{fa} is for a naive factor analytic correlation structure.
#'
#' @param factor Data factor.
#' @param m Number of factors.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
fa <- function(factor, m) {
  levels <- levels(factor)
  p <- length(levels)
  nparams <- p * m
  return(list(base = NA,
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = rep("L", nparams),
              indexCor = 1:nparams,
              part = "var"))
}

#' \code{exfa} is for an extended (sparse) factor analytic correlation structure.
#'
#' @param factor Data factor.
#' @param m Number of factors.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
exfa <- function(factor, m) {
  levels <- levels(factor)
  p <- length(levels)
  nparams <- p * m
  return(list(base = NA,
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = rep("L", nparams),
              indexCor = 1:nparams,
              part = "var"))
}

#' \code{cors} is for a simple correlation structure.
#'
#' @param factor Data factor.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
cors <- function(factor) {
  levels <- levels(factor)
  p <- length(levels)
  nparams <- 1
  return(list(base = NA,
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = rep("R", nparams),
              indexCor = list(1:((p^2 - p) / 2)),
              part = "var"))
}

#' \code{id} is for a diagonal (identity) correlation structure.
#'
#' @param factor Data factor.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
id <- function(factor) {
  levels <- levels(factor)
  p <- length(levels)
  nparams <- 0
  return(list(base = diag(p),
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = NA,
              indexCor = NA,
              part = "nvar"))
}

#' \code{gk} is for a Gaussian kernel correlation structure.
#'
#' @param factor Data factor.
#' @param base Distance matrix.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
gk <- function(factor, base) {
  if (!(all(diag(base) == 0))) {
    stop("Non-zero diagonal entries in distance matrix. Has it been mistakingly inverted?")
  }
  levels <- levels(factor)
  p <- length(levels)
  nparams <- 1
  return(list(base = base[levels, levels],
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = "B",
              indexCor = 0,
              part = "var"))
}

#' \code{lk} is for a linear kernel correlation structure.
#'
#' @param factor Data factor.
#' @param base Correlation matrix or kinship.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
lk <- function(factor, base) {
  levels <- levels(factor)
  p <- length(levels)
  nparams <- 0
  return(list(base = base[levels, levels],
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = NA,
              indexCor = NA,
              part = "nvar"))
}

#' stack
#'
#' @param formula The formula for the stacked correlation matrices
#' @param ... 
#'
#' @returns A structure suitable for the stacked correlation matrices
#' @keywords internal
stack <- function(formula, ...) {
  data <- list(...)$data
  calls <- strsplit(as.character(formula[2]), "\\s*\\+\\s*")[[1]]
  results <- lapply(calls, function(call) {
    eval(parse(text = call), envir = data)
  })
  bases <- lapply(results, `[[`, "base")
  inv <- lapply(bases, function(base) attributes(base)$INVERSE)
  if (any(unlist(inv))) {
    stop("Base matrices for linear kernels that are part of a stack should not be inverted!")
  }
  levels <- unique(lapply(results, `[[`, "levels")) # Every element of this list should be the same, so unique() should return a list of length 1
  stopifnot(length(levels) == 1)
  nlevels <- unique(vapply(results, `[[`, numeric(1), "nlevels")) # Every element of this vector should be the same, so unique() should return a single number
  stopifnot(length(nlevels) == 1)
  nparams <- lapply(results, `[[`, "nparams")
  types <- lapply(results, `[[`, "types")
  indexCor <- lapply(results, `[[`, "indexCor")
  # part <- vapply(results, `[[`, character(1), "part")
  stackResults <- list(base = bases,
                       levels = unlist(levels),
                       nlevels = nlevels,
                       nparams = nparams,
                       types = types,
                       indexCor = indexCor,
                       part = "var", # Stacks are by definition in the variable part of the random effect
                       sizeOfStack = length(calls))
  return(stackResults)
}

#' kmmr
#'
#' @param fixed Formular for the fixed effects.
#' @param random Formula for the random effects.
#' @param residual Formula for the residual effects.
#' @param data Dataframe.
#' @param maxit Maximum number of iterations.
#' @param stepSize Numeric vector indicating the step sizes for parameter updates in each iteration.
#' @param avInfReg Numeric vector indicating the penalty to be applied to the average information matrix in each iteration. Using \code{1.1} inflates the diagonal
#' elements of the AI matrix by 10% before inverting.
#' @param weights Optional weights, i.e., standard errors or numbers of samples.
#' @param init Initial values.
#' @param tolConvLL Tolerance for detecting convergence.
#' @param isFixed Should any of the parameters be fixed?
#' @param verbose Verbosity.
#' @param sparse Optional argument for forcing the sparse or dense versions of the solver.
#'
#' @returns A list with the results.
#' @export
kmmr <- function(fixed = NULL, random = NULL, residual = ~ units, data = NULL,
                 maxit = 20, stepSize = c(0.5, 0.75, rep(1, maxit - 2)), avInfReg = rep(1, maxit),
                 weights = NULL, init = NULL, tolConvLL = 0.001, isFixed = NULL, verbose = TRUE,
                 sparse) {
  
  # Basic argument checking:
  if (is.null(fixed)) {
    stop("At least one fixed effect must be specified.")
  }
  if (is.null(random)) {
    stop("At least one random effect must be specified.")
  }
  if (residual != ~ units) {
    stop("Only the '~ units' residual structure has currently been implemented.")
  }
  if (is.null(data)) {
    stop("The 'data' argument must be specified.")
  }
  if (!(is.numeric(maxit))) {
    stop("The 'maxit' argument must be a positive integer.")
  } else if (length(maxit) > 1) {
    stop("The 'maxit' argument must have length 1.")
  }
  if (!is.numeric(stepSize)) {
    stop("The 'stepSize' argument must be a numeric vector.")
  } else if (length(stepSize) != maxit) {
    stop("The 'stepSize' argument must be of length 'maxit'.")
  }
  if (!(is.character(weights) | is.null(weights))) {
    stop("The 'weights' argument must be NULL or a character.")
  } else if (length(weights) > 1) {
    stop("The 'weights' argument must be a character vector of length 1.")
  }
  if (!(is.null(init) | is.numeric(init))) {
    stop("The 'init' argument must be NULL or numeric.")
  }
  if (!is.numeric(tolConvLL)) {
    stop("The 'tolConvLL' argument must be numeric.")
  } else if (length(tolConvLL) > 1) {
    stop("The 'tolConvLL' argument must have length 1.")
  }
  if (!(is.logical(isFixed) | is.null(isFixed))) {
    stop("The 'isFixed' argument must be NULL or a logical.")
  }
  if (!is.logical(verbose)) {
    stop("The 'verbose' argument must be a logical.")
  } else if (length(verbose) > 1) {
    stop("The 'verbose' argument must have length 1.")
  }
  if (is.numeric(weights)) {
    weightsIsMatrix <- is.matrix(weights)
  } else if (is.character(weights)) {
    weightsIsMatrix <- FALSE
    if (!(weights %in% colnames(data))) {
      stop("The 'weights' argument does not match a column name of the dataframe.")
    }
  } else if (!is.null(weights)) {
    stop("The 'weights' argument must be the name of a column in the dataframe, numeric vector of weights, or numeric matrix of weights.")
  }
  
  # Getting the name of the response column:
  response <- as.character(fixed[2])
  
  # Getting rid of any missing values and adding a unit factor:
  data <- na.omit(data)
  data$Units <- factor(paste0("U", 1:nrow(data)), paste0("U", 1:nrow(data)))
  
  # Fixed effects:
  X <- sparse.model.matrix(revInteractions(fixed), data)
  
  # Random effects:
  res <- strsplit(as.character(random[2]), "\\s*\\+\\s*(?=re\\()", perl = TRUE)[[1]] # split the full random argument into separate re() calls.
  nU <- length(res)
  stepSize <- stepSize^nU
  
  # Lists to store all specifications:
  Z <- indices <- u.applyvars <- u.types <- u.nlevels <- u.forFactor <-
    u.indexCor <- bases <- struc <- u.nparams <- i.var <- i.nvar  <-
    u.isStack <- u.nparams.list <- vector("list", nU)
  params <- vector("list", nU + 1)
  i <- 1
  for (i in 1:nU) {
    expr <- parse(text = res[i])[[1]]
    expr$data <- data
    spec <- eval(expr)
    Z[[i]] <- spec$Z
    indices[[i]] <- spec$indices
    u.applyvars[[i]] <- spec$applyvars
    u.types[[i]] <- spec$types
    u.nlevels[[i]] <- spec$nlevels
    u.forFactor[[i]] <- spec$forFactor
    u.indexCor[[i]] <- spec$indexCor
    bases[[i]] <- spec$bases
    struc[[i]] <- spec$struc
    u.nparams[[i]] <- spec$nparams
    i.var[[i]] <- spec$i.var
    i.nvar[[i]] <- spec$i.nvar
    params[[i]] <- sprintf("re%s(%s)", i, spec$params)
    u.isStack[[i]] <- spec$isStack
    u.nparams.list[[i]] <- spec$nparams.list
  }
  
  # Shifting indices if there are multiple random effects so we can use a single big vector k:
  if (nU > 1) {
    for (i in 2:nU) {
      # indices[[i]] <- lapply(indices[[i]], function(x) x + u.nparams[[i - 1]])
      indices[[i]] <- lapply(indices[[i]], function(x) x + sum(unlist(u.nparams[1:(i - 1)])))
    }
  }
  
  # Residual effects, just assume there is one simple iid effect for now:
  nR <- 1
  S <- vector("list", nR)
  S[[1]] <- sparse.model.matrix(~ -1 + Units, data)
  ur.types <- c(unlist(u.types), "V") # This assumes just a single residual variance.
  params[[nU + 1]] <- "residual(vars{units})"
  
  # Response matrix:
  y <- Matrix(data[[response]], sparse = FALSE, doDiag = FALSE)
  
  if (is.null(weights)) {
    useWeights <- FALSE
    weights <- NA
  } else {
    useWeights <- TRUE
    if (weightsIsMatrix) {
      weights <- chol(weights) # Let's not worry about the full covariance matrix yet, in this case 'weights' would be omega^-1
    } else {
      # For the diagonal approximation, 'weights' is a vector containing the diagonal elements of omega^-1, see Damesa et al. 2017, https://doi.org/10.2134/agronj2016.07.0395
      # In this case we can estimate the inverse residual correlation matrix Ri if needed, or just use an identity matrix for independent residuals so that Ri = I.
      # Then multiply on the left and right with a diagonal matrix of the square roots of the original weights: weights %*% Ri %*% weights.
      # Or maybe alternatively, solve(outer(diag(weights), diag(weights)) * Ri) to obtain the residual covariance matrix.
      # (outer(diag(weights), diag(weights)) * Ri is the inverse residual covariance matrix).
      # Note that any residual variances must be fixed to 1.0
      if (is.numeric(weights)) {
        weights <- diag(sqrt(weights))
      } else {
        weights <- Matrix(diag(sqrt(data[[weights]])))
      }
    }
  }
  
  # Initial values:
  if (is.null(init)) {
    init <- rep(0.3, length(ur.types))
  } else if (length(init) != length(ur.types)) {
    stop("Length of 'init' does not match the number of parameters.")
  } else if (any(is.na(init))) {
    init[is.na(init)] <- 0.3
  }
  
  # Fixed values:
  if (is.null(isFixed)) {
    isFixed <- rep(0, length(init))
  } else if (length(isFixed) != length(ur.types)) {
    stop("Length of 'isFixed' does not match the number of parameters.")
  } else {
    isFixed <- as.numeric(isFixed)
  }
  
  # Fitting:
  # browser()
  sparse = FALSE
  result <- solveLMM_R(y = y, X = X, Z = Z, S = S, init = init, weights = weights, useWeights = useWeights,
                       maxit = maxit, stepSize = stepSize, avInfReg = avInfReg, verbose = verbose,
                       tolConvLL = tolConvLL, indices = indices, u.applyvars = u.applyvars,
                       ur.types = ur.types, u.nlevels = u.nlevels, u.nparams = u.nparams,
                       u.forFactor = u.forFactor, u.indexCor = u.indexCor,
                       bases = bases, struc = struc, i.var = i.var, i.nvar = i.nvar,
                       isFixed = isFixed, u.isStack = u.isStack, u.nparams.list = u.nparams.list, sparse = sparse)
  result$params <- unlist(params)
  class(result) <- "kmmrFit"
  return(result)
}


