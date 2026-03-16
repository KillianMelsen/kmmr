#' revInteractions
#'
#' @param form Formula.
#'
#' @returns The formula with kronecker products/interaction terms reversed.
#' @keywords internal
#' @importFrom stats terms as.formula
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
#' @importFrom stats na.omit
#' @importFrom Matrix sparse.model.matrix
re <- function(formula, vars, ...) {
  data <- list(...)$data
  calls <- strsplit(as.character(formula[2]), ":")[[1]] # vector of function calls for each structure/factor.
  results <- lapply(calls, function(call) {
    eval(parse(text = call), envir = data)
  })
  
  struc <- gsub("(.*)\\(.*$", "\\1", calls) # What are the correlation structures for this factor?
  factorNames <- gsub(".*\\(([^,\\)]+).*$", "\\1", calls) # Which factors are we evaluating?
  
  bases <- lapply(results, `[[`, "base")
  levels <- lapply(results, `[[`, "levels")
  nlevels <- vapply(results, `[[`, numeric(1), "nlevels")
  nparams <- vapply(results, `[[`, numeric(1), "nparams")
  part <- vapply(results, `[[`, character(1), "part")
  types <- unlist(lapply(results, `[[`, "types"))
  indexCor <- unlist(lapply(results, `[[`, "indexCor"), recursive = FALSE)
  forFactor <- unlist(lapply(seq_along(results), function(i) {rep(i, results[[i]]$nparams)}))
  
  if (as.character(vars[2]) == "1") {
    applyvars <- 0
    nvars <- 1
  } else {
    applyvars <- which(factorNames %in% strsplit(as.character(vars[2]), "[:]")[[1]])
    nvars <- prod(nlevels[applyvars])
  }
  
  # Determine which part of the random effect is variable and which part is non-variable:
  i.var <- setdiff(unique(c(applyvars, which(part == "var"))), 0)
  i.nvar <- setdiff(1:length(calls), i.var)
  if (!(all(i.var == 1:length(i.var)))) {
    stop(paste0("In '~ re(~ ", as.character(formula[2]), ", vars = ~ ", as.character(vars[2]), ")', linear kernels to which no variances are applied must be placed in the right part of the kronecker product."))
  }
  
  # Now we can add the variances to all vectors:
  types <- c(rep("V", nvars), na.omit(types))
  forFactor <- c(rep(0, nvars), forFactor)
  nparams <- c(nvars, nparams)
  indexCor <- c(rep(0, nvars), indexCor[!is.na(indexCor)])
  
  indices <- vector("list", 1 + length(calls)) # We will fill this later after we figure out all factor structures.
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
  
  return(list(Z = Z,
              indices = indices,
              applyvars = applyvars,
              types = types,
              nlevels = nlevels,
              forFactor = forFactor,
              indexCor = indexCor,
              bases = bases,
              struc = struc,
              nparams = sum(nparams),
              i.var = i.var,
              i.nvar = i.nvar))
}

#' corg
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

#' cors
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

#' id
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

#' gk
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

#' lk
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

#' kmmr
#'
#' @param fixed Formular for the fixed effects.
#' @param random Formula for the random effects.
#' @param residual Formula for the residual effects.
#' @param data Dataframe.
#' @param maxit Maximum number of iterations.
#' @param stepsize Stepsizes for parameters updates.
#' @param weights Optional weights, i.e., standard errors or numbers of samples.
#' @param init Initial values.
#' @param tolConvLL Tolerance for detecting convergence.
#' @param isFixed Should any of the parameters be fixed?
#' @param verbose Verbosity.
#'
#' @returns A list with the results.
#' @export
#' @importFrom stats na.omit
#' @importFrom Matrix sparse.model.matrix
kmmr <- function(fixed = NULL, random = NULL, residual = ~ units, data = NULL,
                 maxit = 20, stepsize = c(0.5, 0.75, rep(1, maxit - 2)), weights = NULL,
                 init = NULL, tolConvLL = 0.001, isFixed = NULL, verbose = TRUE) {
  
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
  if (!is.numeric(stepsize)) {
    stop("The 'stepsize' argument must be a numeric vector.")
  } else if (length(stepsize) != maxit) {
    stop("The 'stepsize' argument must be of length 'maxit'.")
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
  
  # Getting the name of the response column:
  response <- as.character(fixed[2])
  
  # Getting rid of any missing values and adding a unit factor:
  data <- na.omit(data)
  data$Units <- factor(paste0("U", 1:nrow(data)), paste0("U", 1:nrow(data)))
  
  # Fixed effects:
  X <- sparse.model.matrix(revInteractions(fixed), data)
  
  # Random effects:
  res <- strsplit(as.character(random[2]), "[+]")[[1]]
  nU <- length(res)
  
  # Lists to store all specifications:
  Z <- indices <- u.applyvars <- u.types <- u.nlevels <- u.forFactor <- u.indexCor <- bases <- struc <- u.nparams <- i.var <- i.nvar <- vector("list", nU)
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
  }
  
  # Shifting indices if there are multiple random effects so we can use a single big vector k:
  if (nU > 1) {
    for (i in 2:nU) {
      indices[[i]] <- lapply(indices[[i]], function(x) x + u.nparams[[i - 1]])
    }
  }
  
  # Residual effects, just assume there is one simple iid effect for now:
  nR <- 1
  S <- vector("list", nR)
  S[[1]] <- sparse.model.matrix(~ -1 + Units, data)
  ur.types <- c(unlist(u.types), "V") # This assumes just a single residual variance.
  
  # Response matrix:
  y <- matrix(data[[response]])
  
  if (is.null(weights)) {
    useH <- FALSE
    H <- NA
  } else {
    useH <- TRUE
    stop("Weights are not yet implemented.")
    # Do something here.
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
  result <- solveLMM_R(y = y, X = X, Z = Z, S = S, init = init, H = H, useH = useH,
                       maxit = maxit, stepsize = stepsize, verbose = verbose,
                       tolConvLL = tolConvLL, indices = indices, u.applyvars = u.applyvars,
                       ur.types = ur.types, u.nlevels = u.nlevels, u.forFactor = u.forFactor,
                       u.indexCor = u.indexCor, bases = bases, struc = struc, i.var = i.var,
                       i.nvar = i.nvar, isFixed = isFixed)
  return(result)
}


