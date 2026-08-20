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
  # rotmat <- Reduce("%*%", givens)
  return(L.rotated)
}

getInit_fa <- function(factor, m, data, response) {
  idCols <- setdiff(colnames(data), c(factor, response, "units"))
  dataWide <- reshape(data[, setdiff(colnames(data), "units")], direction = "wide", idvar = idCols, timevar = factor, v.names = response)
  dataWide <- dataWide[, setdiff(colnames(dataWide), idCols)]
  colnames(dataWide) <- gsub(paste0("^", response, "\\."), "", colnames(dataWide))
  dataWide <- as.matrix(dataWide)[, levels(data[[factor]])]
  covmat <- cov(dataWide, use = "pairwise.complete.obs") # Normal covariance matrix if there are no missing values (equal to use = "complete").
  if (any(is.na(covmat))) { # If we have pairs with no overlap and thus NAs, we need to impute using column means.
    means <- colMeans(dataWide, na.rm = TRUE)
    naCoords <- which(is.na(dataWide), arr.ind = TRUE)
    dataWide[naCoords] <- means[naCoords[, "col"]]
    covmat <- cov(dataWide, use = "pairwise.complete.obs")
    if (min(eigen(covmat, only.values = TRUE)$values) <= 0) { # Still check for PDness
      covmat <- as.matrix(nearPD(covmat, tol)$mat)
    }
  } else if (min(eigen(covmat, only.values = TRUE)$values) <= 0) { # If we don't have any NAs, we still need to check whether covmat is PD.
    covmat <- as.matrix(nearPD(covmat, tol)$mat)
  }
  fit.factanal <- factanal(covmat = covmat, factors = m)
  L <- fit.factanal$loadings[,, drop = FALSE]
  if (ncol(L) > 1) {
    L <- rotateInitL(L)
  }
  L <- as.numeric(L)
  return(list(init = L, isFixed = L == 0))
}

getInit_exfa <- function(factor, m, data, response) {
  return(getInit_fa(factor, m, data, response))
}

getInit_gk <- function(factor, base, data, response) {
  return(list(init = 0.15, isFixed = FALSE))
}

getInit_cors <- function(factor, data, response) {
  return(list(init = 0.5, isFixed = FALSE))
}

getInit_corg <- function(factor, data, response) {
  p <- length(levels(data[[factor]]))
  ncors <- (p^2 - p) / 2
  return(list(init = rep(0.5, ncors), isFixed = rep(FALSE, ncors)))
}

getInit_ar1 <- function(factor, data, response) {
  return(list(init = 0.1, isFixed = FALSE))
}

getInit_stack <- function(formula, data, response) {
  callsInStack <- strsplit(as.character(formula[2]), "\\s*\\+\\s*")[[1]]
  resultStack <- lapply(callsInStack, function(call) {
    if (!(gsub("(.*)\\(.*", "\\1", call) %in% c("lk", "id"))) { # We don't need to do anything if the structure is id or lk as these have no parameters.
      call <- gsub("(.*\\()([^,]*)(,?.*)\\)", "getInit_\\1'\\2'\\3, data = data, response = response)", call)
      result <- eval(parse(text = call)[[1]])
    }
  })
  initStack <- unlist(lapply(resultStack, `[[`, "init"))
  isFixedStack <- unlist(lapply(resultStack, `[[`, "isFixed"))
  return(list(init = initStack, isFixed = isFixedStack))
}

#' getInit
#'
#' @param response Name of the response variable.
#' @param call \code{re()} call as a string.
#' @param data Input dataframe.
#' @param init Current initial values vector.
#' @param isFixed Current vector for fixing parameters.
#' @param indices Parameter indices for this random effect.
#'
#' @returns Modified init and isFixed vectors.
#' @keywords internal
getInit <- function(response, call, data, init, isFixed, indices) {
  calls <- strsplit(gsub("re\\(\\s*~\\s*(.*)\\s*,.*", "\\1", call), ":")[[1]]
  for (j in 1:length(calls)) {
    currentCall <- calls[j]
    if (grepl("stack\\(.*", currentCall)) { # If it's a stacked structure.
      currentCall <- gsub("(.*\\()(.*)\\)", "getInit_\\1\\2, data = data, response = response)", currentCall)
      result <- eval(parse(text = currentCall)[[1]])
    } else if (!(gsub("(.*)\\(.*", "\\1", currentCall) %in% c("lk", "id", "exfalk"))) { # We don't need to do anything if the structure is id, lk, or exfalk as these have no parameters.
      currentCall <- gsub("(.*\\()([^,]*)(,?.*)\\)", "getInit_\\1'\\2'\\3, data = data, response = response)", currentCall)
      result <- eval(parse(text = currentCall)[[1]])
      
      # 1 + j so we skip indices for the variances.
      # We leave init values that are not NA. These have been supplied by the user. We will only
      # set initial values that were not yet set by the user.
      init[indices[[1 + j]]][is.na(init[indices[[1 + j]]])] <- result$init[is.na(init[indices[[1 + j]]])]
      isFixed[indices[[1 + j]]][is.na(isFixed[indices[[1 + j]]])] <- result$isFixed[is.na(isFixed[indices[[1 + j]]])]
    }
  }
  return(list(init = init, isFixed = isFixed))
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
    stackInCalls <- grepl("stack\\(.*", call)
    if (stackInCalls) {
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
  # If we have an exfa structure we need to add some levels for the factors:
  if (any(struc == "exfa")) {
    levelsExtended <- levels # Make a copy of the list of levels.
    for (f in which(struc == "exfa")) {
      m <- nparams[1 + f] / length(levels[[f]]) # Add 1 because the first entry is the number of variances.
      levelsExtended[[f]] <- c(levelsExtended[[f]], sprintf("%s-F%d", factorNames[f], 1:m))
    }
    superlevels <- do.call(paste, c(rev(expand.grid(rev(levelsExtended))), sep = ":"))
  } else {
    superlevels <- do.call(paste, c(rev(expand.grid(rev(levels))), sep = ":"))
  }
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
    if (nparams[i + 1] > 0) {
      if (isStack[i]) {
        for (s in 1:sos) {
          p <- nparams.list[[i + 1]][s]
          params[h:(h - 1 + p)] <- addParamLabels(structure = struc[[i]][s], levels = levels[[i]], factorName = factorNames[i], p = p, s = s)
          h <- h + p
        }
      } else {
        p <- nparams[i + 1]
        params[h:(h - 1 + p)] <- addParamLabels(structure = struc[i], levels = levels[[i]], factorName = factorNames[i], p = p)
        h <- h + p
      }
    }
  }
  # This should probably be wrapped in a function...
  # for (i in 1:length(factorNames)) {
  #   if (nparams[i + 1] > 0) { # + 1 because we can skip the variances
  #     if (isStack[i]) {
  #       for (s in 1:sos) {
  #         p <- nparams.list[[i + 1]][s]
  #         structure <- struc[[i]][s]
  #         if (structure == "corg") {
  #           labels <- t(outer(levels[[i]], levels[[i]], function(x, y) paste(x, y, sep = ":")))[lower.tri(matrix(1, nlevels[i], nlevels[i]))]
  #           params[h:(h - 1 + p)] <- sprintf("corg{%s-%s|stack-%d}", factorNames[i], labels, s)
  #           h <- h + p
  #         } else if (structure == "cors") {
  #           params[h] <- sprintf("cors{%s|stack-%d}", factorNames[i], s)
  #           h <- h + 1
  #         } else if (structure == "gk") {
  #           params[h] <- sprintf("bw{%s|stack-%d}", factorNames[i], s)
  #           h <- h + 1
  #         } else if (structure %in% c("fa", "exfa")) {
  #           m <- p / nlevels[i]
  #           labels <- paste(rep(levels[[i]], m), rep(sprintf("F%d", 1:m), each = nlevels[i]), sep = ":")
  #           params[h:(h - 1 + p)] <- sprintf("loading{%s-%s|stack-%d}", factorNames[i], labels, s)
  #           h <- h + p
  #         } else if (structure == "ar1") {
  #           params[h] <- sprintf("ar1{%s}", factorNames[i])
  #           h <- h + 1
  #         }
  #       }
  #     } else {
  #       p <- nparams[i + 1]
  #       structure <- struc[i]
  #       if (structure == "corg") {
  #         labels <- t(outer(levels[[i]], levels[[i]], function(x, y) paste(x, y, sep = ":")))[lower.tri(matrix(1, nlevels[i], nlevels[i]))]
  #         params[h:(h - 1 + p)] <- sprintf("corg{%s-%s}", factorNames[i], labels)
  #         h <- h + p
  #       } else if (structure == "cors") {
  #         params[h] <- sprintf("cors{%s}", factorNames[i])
  #         h <- h + 1
  #       } else if (structure == "gk") {
  #         params[h] <- sprintf("bw{%s}", factorNames[i])
  #         h <- h + 1
  #       } else if (structure %in% c("fa", "exfa")) {
  #         m <- p / nlevels[i]
  #         labels <- paste(rep(levels[[i]], m), rep(sprintf("F%d", 1:m), each = nlevels[i]), sep = ":")
  #         params[h:(h - 1 + p)] <- sprintf("loading{%s-%s}", factorNames[i], labels)
  #         h <- h + p
  #       } else if (structure == "ar1") {
  #         params[h] <- sprintf("ar1{%s}", factorNames[i])
  #         h <- h + 1
  #       }
  #     }
  #   }
  # }
  
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
              params = params,
              factorNames = factorNames))
}

addParamLabels <- function(structure, levels, factorName, p, s = NULL) {
  # Making the standard labels:
  if (structure == "corg") {
    labels <- t(outer(levels, levels, function(x, y) paste(x, y, sep = ":")))[lower.tri(matrix(1, length(levels), length(levels)))]
    labels <- sprintf("corg{%s-%s}", factorName, labels)
  } else if (structure == "cors") {
    labels <- sprintf("cors{%s}", factorName)
  } else if (structure == "gk") {
    labels <- sprintf("bw{%s}", factorName)
  } else if (structure %in% c("fa", "exfa")) {
    m <- p / length(levels)
    labels <- paste(rep(levels, m), rep(sprintf("F%d", 1:m), each = length(levels)), sep = ":")
    labels <- sprintf("loading{%s-%s}", factorName, labels)
  } else if (structure == "ar1") {
    labels <- sprintf("ar1{%s}", factorName)
  }
  
  # Adding the stack label if required:
  if (!is.null(s)) {
    labels <- gsub("(.*)}", sprintf("\\1|stack-%d}", s), labels)
  }
  
  return(labels)
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
  return(list(base = Diagonal(p),
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

#' \code{exfalk} is for an extended factor analytic approximation of a linear kernel.
#' Only used for approximating a fully dense kinship matrix.
#'
#' @param factor Data factor.
#' @param setup Setup returned from the \code{exfalkSetup()} function.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
exfalk <- function(factor, setup) {
  levels <- c(levels(factor), sprintf("%s-F%d", deparse(substitute(factor)), 1:setup$m))
  p <- length(levels)
  nparams <- 0
  PSIi_L <- setup$PSIi %*% setup$L
  base <- forceSymmetric(rbind(cbind(setup$PSIi, PSIi_L),
                               cbind(t(PSIi_L), Diagonal(setup$m) + crossprod(setup$L, PSIi_L))))
  rownames(base) <- colnames(base) <- c(rownames(setup$L), sprintf("%s-F%d", deparse(substitute(factor)), 1:setup$m))
  if (!all.equal(rownames(base), levels)) {
    base <- base[levels, levels]
  }
  return(list(base = base,
              levels = levels,
              nlevels = p,
              nparams = nparams,
              types = NA,
              indexCor = NA,
              part = "nvar"))
}

#' \code{ar1} is for an order 1 autoregressive correlation structure
#'
#' @param factor Data factor.
#'
#' @returns List with all information on this correlation structure.
#' @keywords internal
ar1 <- function(factor) {
  levels <- levels(factor)
  indices <- as.numeric(levels)
  D <- Matrix(abs(outer(indices, indices, FUN = "-")), dimnames = list(levels, levels))
  return(list(base = D,
              levels = levels,
              nlevels = length(levels),
              nparams = 1,
              types = "AR",
              indexCor = NA,
              part = "var"))
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

atEach <- function(factor, formula, levels = NULL, ...) {
  data <- list(...)$data
  if (is.null(levels)) {
    levels <- levels(data[[factor]])
  } else if (any(!(levels %in% levels(data[[factor]])))) {
    stop("The 'levels' argument of one of the atEach() calls contains levels not present in the data...")
  }
  specs <- vector("list", length(levels))
  l <- 1
  for (lvl in levels) {
    rowIndices <- which(data[[factor]] == lvl)
    dataSubset <- droplevels(data[rowIndices,])
    expr <- parse(text = as.character(formula[2]))[[1]]
    expr$data <- dataSubset
    specs[[l]] <- eval(expr)
    expandedZ <- Matrix(0, nrow(data), ncol(specs[[l]]$Z), sparse = TRUE, doDiag = FALSE)
    expandedZ[rowIndices,] <- specs[[l]]$Z
    colnames(expandedZ) <- colnames(specs[[l]]$Z)
    specs[[l]]$Z <- expandedZ
    names(specs)[l] <- sprintf("%s-%s", factor, lvl)
    l <- l + 1
  }
  return(specs)
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
#' @param timesToBoundary How many times is a parameter or set of loadings allowed to cross the boundary before it is fixed?
#' @param getStartingValues Should the model be fitted, or should suggested starting values be returned?
#' @param super Should a supernodal (\code{super = TRUE}) or simplicial Cholesky factor be computed?
#'
#' @returns A list with the results.
#' @export
kmmr <- function(fixed = NULL, random = NULL, residual = ~ units, data = NULL,
                 maxit = 20, stepSize = NULL, avInfReg = NULL,
                 weights = NULL, init = NULL, tolConvLL = 0.001, isFixed = NULL, verbose = TRUE,
                 sparse = FALSE, timesToBoundary = 5, getStartingValues = FALSE, super = TRUE) {
  
  # Basic argument checking:
  if (is.null(fixed)) {
    stop("At least one fixed effect must be specified.")
  }
  if (is.null(random)) {
    stop("At least one random effect must be specified.")
  }
  # if (residual != ~ units) {
  #   stop("Only the '~ units' residual structure has currently been implemented.")
  # }
  if (is.null(data)) {
    stop("The 'data' argument must be specified.")
  }
  if (!(is.numeric(maxit))) {
    stop("The 'maxit' argument must be a positive integer.")
  } else if (length(maxit) > 1) {
    stop("The 'maxit' argument must have length 1.")
  }
  if (!(is.numeric(stepSize) | is.null(stepSize))) {
    stop("The 'stepSize' argument must be a numeric vector or NULL.")
  } else if (is.numeric(stepSize) & length(stepSize) > maxit) {
    stop("A numeric 'stepSize' argument must be of length 'maxit' or smaller.")
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
  
  # Getting rid of any missing values:
  data <- na.omit(data)
  
  # Fixed effects:
  X <- sparse.model.matrix(revInteractions(fixed), data)
  
  # Random effects and residual effect:
  # Short-hand notation for the residuals: 'residual = ~ units'
  # This results in 'residual = ~ re(~ id(units), vars = ~ 1)' where
  # units is an internally made factor in the dataframe where each row corresponds
  # to one level. So it results in iid residuals.
  if (as.character(residual[2]) == "units") { # Short-hand 'residual = ~ units'...
    residual <- ~ re(~ id(units), vars = ~ 1) #...becomes 'residual = ~ re(~ id(units), vars = ~ 1)', that is, iid residuals.
    data$units <- factor(sprintf("U%d", 1:nrow(data)), sprintf("U%d", 1:nrow(data)), sprintf("U%d", 1:nrow(data)))
  }
  
  # S <- Diagonal(n.units)
  res <- strsplit(as.character(random[2]), "\\s*\\+\\s*(?=re\\()", perl = TRUE)[[1]] # split the full random argument into separate re() calls.
  nU <- length(res) # Number of random effects.
  res <- c(res, as.character(residual[2])) # Residual effect
  n.atEach <- sum(grepl("^atEach.*", res))
  
  ## Working on atEach() ====
  if (n.atEach > 0) {
    specs.atEach <- vector("list", n.atEach)
    for (i in 1:n.atEach) {
      expr <- parse(text = res[which(grepl("^atEach.*", res))[i]])[[1]]
      expr$data <- data
      specs.atEach[[i]] <- eval(expr)
    }
    # We need to calculate how many random effects we now actually have (+1 for the residual effect):
    nU <- 0
    i.atEach <- 1 # Index to keep track of which atEach() call we are.
    for (i in 1:(length(res) - 1)) {
      if (grepl("^atEach.*", res[i])) {
        nU <- nU + length(specs.atEach[[i.atEach]])
        i.atEach <- i.atEach + 1
      } else {
        nU <- nU + 1
      }
    }
  }
  
  # Lists to store all specifications, note that the very last entry is reserved for the residual effect:
  # The length of these lists already takes additional random effects due to atEach() calls into account.
  Z <- indices <- u.applyvars <- u.types <- u.nlevels <- u.forFactor <-
    u.indexCor <- bases <- struc <- u.nparams <- i.var <- i.nvar  <-
    u.isStack <- u.nparams.list <- params <- factorNames <- vector("list", nU + 1)
  j <- 1 # Index to keep track of the actual random effect spec we are working on.
  k <- 1 # Index to keep track of which atEach() spec we are at.
  # for (i in 1:(nU + 1)) {
  for (i in 1:length(res)) { # We now need to iterate over the calls in res.
    if (grepl("^atEach.*", res[i])) { # If we are dealing with an atEach() call that expands into several random effects.
      for (l in 1:length(specs.atEach[[k]])) {
        Z[[j]] <- specs.atEach[[k]][[l]]$Z
        indices[[j]] <- specs.atEach[[k]][[l]]$indices
        u.applyvars[[j]] <- specs.atEach[[k]][[l]]$applyvars
        u.types[[j]] <- specs.atEach[[k]][[l]]$types
        u.nlevels[[j]] <- specs.atEach[[k]][[l]]$nlevels
        u.forFactor[[j]] <- specs.atEach[[k]][[l]]$forFactor
        u.indexCor[[j]] <- specs.atEach[[k]][[l]]$indexCor
        bases[[j]] <- specs.atEach[[k]][[l]]$bases
        struc[[j]] <- specs.atEach[[k]][[l]]$struc
        u.nparams[[j]] <- specs.atEach[[k]][[l]]$nparams
        i.var[[j]] <- specs.atEach[[k]][[l]]$i.var
        i.nvar[[j]] <- specs.atEach[[k]][[l]]$i.nvar
        params[[j]] <- sprintf("atEach%s(re%s(%s|%s))", k, j, specs.atEach[[k]][[l]]$params, names(specs.atEach[[k]])[l])
        u.isStack[[j]] <- specs.atEach[[k]][[l]]$isStack
        u.nparams.list[[j]] <- specs.atEach[[k]][[l]]$nparams.list
        factorNames[[j]] <- specs.atEach[[k]][[l]]$factorNames
        j <- j + 1
      }
      k <- k + 1
    } else { # If we are dealing with an re() call that results in a specification for a single random effect.
      expr <- parse(text = res[i])[[1]]
      expr$data <- data
      spec <- eval(expr)
      Z[[j]] <- spec$Z
      indices[[j]] <- spec$indices
      u.applyvars[[j]] <- spec$applyvars
      u.types[[j]] <- spec$types
      u.nlevels[[j]] <- spec$nlevels
      u.forFactor[[j]] <- spec$forFactor
      u.indexCor[[j]] <- spec$indexCor
      bases[[j]] <- spec$bases
      struc[[j]] <- spec$struc
      u.nparams[[j]] <- spec$nparams
      i.var[[j]] <- spec$i.var
      i.nvar[[j]] <- spec$i.nvar
      if (j <= nU) { # Labels for random effects.
        params[[j]] <- sprintf("re%s(%s)", j, spec$params)
      } else { # Labels for the residual effect.
        params[[j]] <- sprintf("residual(%s)", spec$params)
      }
      u.isStack[[j]] <- spec$isStack
      u.nparams.list[[j]] <- spec$nparams.list
      factorNames[[j]] <- spec$factorNames
      j <- j + 1
    }
  }
  
  # Shifting indices so we can use a single big vector k:
  for (i in 2:(nU + 1)) {
    indices[[i]] <- lapply(indices[[i]], function(x) x + sum(unlist(u.nparams[1:(i - 1)])))
  }
  
  # Vector of parameter types
  ur.types <- unlist(u.types)
  
  # Response matrix:
  y <- Matrix(data[[response]], sparse = FALSE, doDiag = FALSE)
  
  if (is.null(weights)) {
    useWeights <- FALSE
    weights <- NA
  } else {
    useWeights <- TRUE
    # For the diagonal approximation, 'weights' is the name of the column of the dataframe containing the diagonal elements of omega^-1, see Damesa et al. 2017, https://doi.org/10.2134/agronj2016.07.0395
    # In this case we can estimate the inverse residual correlation matrix Ri if needed, or just use an identity matrix for independent residuals so that Ri = I.
    # Then multiply on the left and right with a diagonal matrix of the square roots of the original weights: weights %*% Ri %*% weights.
    # Note that any residual variances must be fixed to 1.0
    weights <- Diagonal(nrow(data), sqrt(data[[weights]])) # Will be a ddiMatrix
  }
  
  # Initial values:
  factorNames <- unique(unlist(factorNames))
  data <- data[, c(factorNames, response)] # Dropping any irrelevant columns to avoid issues with getting initial values.
  if (is.null(init) | any(is.na(init))) {
    if (is.null(init)) {
      init <- rep(NA, length(ur.types))
    }
    if (is.null(isFixed)) {
      isFixed <- rep(NA, length(ur.types))
    }
    for (i in 1:(nU + 1)) {
      # Get initial variances:
      formula <- as.formula(paste0(response, gsub(".*vars\\s*=\\s*(~.*)\\)", "\\1", res[i])))
      if (formula == as.formula(paste0(response, "~ 1"))) {
        vars <- var(data[[response]])
      } else {
        # vars <- diag(cov(unstack(data, formula), use = "complete"))
        vars <- aggregate(data, formula, var)[[response]]
      }
      init[indices[[i]][[1]]] <- vars / 20
      isFixed[indices[[i]][[1]]] <- FALSE
      
      # Now do the initial values for the remaining structures in this random effect (or residual effect if 1 == nU + 1):
      ret <- getInit(response = response, call = res[i], data = data, init = init, isFixed = isFixed, indices = indices[[i]])
      init <- ret$init
      isFixed <- ret$isFixed
    }
  }
  # Setting the AI matrix penalties if not set by the user:
  # Also adding some penalty to the average information matrix if factor analytic structures are present:
  if (is.null(avInfReg)) {
    if (any(unlist(struc) %in% c("fa", "exfa"))) {
      avInfReg <- c(1.1, 1.05, rep(1, maxit - 2))
    } else {
      avInfReg <- rep(1, maxit)
    }
  }
  
  # Setting the step sizes if they were not set by the user:
  if (is.null(stepSize)) {
    stepSize <- c(0.25, 0.5, 0.75, rep(1, maxit - 3))
    stepSize <- stepSize^nU
  } else if (length(stepSize) < maxit) {
    stepSize <- c(stepSize, rep(1, maxit - length(stepSize)))
  }
  
  if (getStartingValues) {
    result <- data.frame(Parameter = unlist(params),
                         Init = init,
                         Fixed = isFixed)
  } else {
    # Fitting:
    # browser()
    result <- solveLMM(y = y, X = X, Z = Z, init = init, weights = weights, useWeights = useWeights,
                       maxit = maxit, stepSize = stepSize, avInfReg = avInfReg, verbose = verbose,
                       tolConvLL = tolConvLL, indices = indices, u.applyvars = u.applyvars,
                       ur.types = ur.types, u.nlevels = u.nlevels, u.nparams = u.nparams,
                       u.forFactor = u.forFactor, u.indexCor = u.indexCor,
                       bases = bases, struc = struc, i.var = i.var, i.nvar = i.nvar,
                       isFixed = isFixed, u.isStack = u.isStack, u.nparams.list = u.nparams.list, sparse = sparse,
                       timesToBoundary = timesToBoundary, super = super)
    result$params <- unlist(params)
    result$trace <- rbind(result$LL[1:ncol(result$k)],
                          stepSize[1:ncol(result$k)],
                          avInfReg[1:ncol(result$k)],
                          result$k)
    rownames(result$trace) <- c("Log-likelihood", "Stepsize", "AI penalty", result$params)
    rownames(result$u) <- unlist(lapply(Z[-length(Z)], colnames)) # Last design matrix is for the residual effect.
    class(result) <- "kmmrFit"
  }
  return(result)
}


