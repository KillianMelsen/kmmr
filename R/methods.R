#' summary.kmmrFit
#'
#' @param object kmmrFit object.
#' @param round Number of decimals to be printed for parameter estimates and SEs.
#' @param ... Extra arguments.
#' 
#' @keywords internal
#' @export
summary.kmmrFit <- function(object, round = 5, ...) {
  estimates <- object$k[, ncol(object$k)]
  ses <- rep(NA, length(estimates))
  ses[which(object$isFixed == 0)] <- sqrt(diag(solve(object$AI[which(object$isFixed == 0), which(object$isFixed == 0)])))
  fixed <- ifelse(object$isFixed == 0, "No", "Yes")
  
  outputData <- list(logLik = object$LL[length(object$LL)],
                     varcomp = data.frame(Parameter = object$params,
                                          Estimate = estimates,
                                          SE = ses,
                                          Fixed = fixed))
  attr(outputData, "round") <- round
  class(outputData) <- c("summary.kmmrFit", "list")
  return(outputData)
  # output <- apply(outputData, 1, function(row) {
  #   sprintf("%-20s %-12.8f %-12.8f %-8s",
  #           row["parameter"],
  #           as.numeric(row["estimate"]),
  #           as.numeric(row["SE"]),
  #           row["isFixed"])
  # })
  # 
  # header <- sprintf("%-20s %-12s %-12s %-8s", "Parameter", "Estimate", "SE", "Fixed?")
  # ansi <- interactive() && !is.null(Sys.getenv("TERM"))
  # if (ansi) {
  #   cat("\033[1m", header, "\033[0m", "\n", sep = "")
  # } else {
  #   cat(header, "\n")
  # }
  # cat(paste(output, collapse = "\n"))
}

#' print.summary.kmmrFit
#'
#' @param x summary.kmmrFit object.
#' @param ... Extra arguments.
#' 
#' @keywords internal
#' @export
print.summary.kmmrFit <- function(x, ...) {
  # Is bold text possible?
  ansi <- interactive() && !is.null(Sys.getenv("TERM"))
  round <- attr(x, "round")
  
  # Formatting parameter estimates, SEs, etc.
  max_len <- max(nchar(x$varcomp$Parameter))
  output <- apply(x$varcomp, 1, function(row) {
    sprintf(gettextf("%%-%ds %%-12.%df %%-12.%df %%-8s", max_len + 5, round, round),
            row["Parameter"],
            as.numeric(row["Estimate"]),
            as.numeric(row["SE"]),
            row["Fixed"])
  })
  
  # Formatting parameter table header:
  header <- sprintf(gettextf("%%-%ds %%-12s %%-12s %%-8s", max_len + 5), "Parameter", "Estimate", "SE", "Fixed?")

  # Printing:
  if (ansi) {
    cat("\033[1m", "Log-likelihood: ", "\033[0m", sprintf("%.3f", x$logLik), "\n\n", sep = "")
    cat("\033[1m", header, "\033[0m", "\n", sep = "")
  } else {
    cat("Log-likelihood: ", sprintf("%.3f", x$logLik), "\n\n", sep = "")
    cat(header, "\n")
  }
  cat(paste(output, collapse = "\n"), "\n")
}




