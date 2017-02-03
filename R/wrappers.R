#' Lightweight model wrappers
#' 
#' @details These functions provide family-specific wrappers for
#'   \code{\link{stan_glm}} and \code{\link{stan_glmer}}. They accept fewer
#'   arguments than \code{stan_glm} and \code{stan_glmer} (e.g., no
#'   \code{subset}, \code{na.action}, \code{constrasts}, etc.), and so in some
#'   cases it is still necessary to use \code{stan_glm} or \code{stan_glmer}
#'   directly. However, these wrappers are nice for three primary reasons.
#'   First, the name of the function is more informative about the model that is
#'   fit. Second, the default priors associated with these wrappers are set on a
#'   family-specific basis (as opposed to \code{stan_glm} and \code{stan_glmer},
#'   which use the same default prior families for all models). Third, these
#'   functions accept both \code{glm}- and \code{glmer}-style formulas.
#' 
#' @name model-wrappers
#' @inheritParams stan_glm
#' @param link The link function to use.
#' @param prior_covariance Ignored if the model has no group-specific terms, 
#'   otherwise passed to \code{\link{stan_glmer}}.
#' 
#' @return The same \link[=stanreg-objects]{stanreg} object that would be 
#'   returned if the model were fit using \code{stan_glm} or \code{stan_glmer}
#'   directly instead of the wrapper.
#'  
NULL

#' @rdname model-wrappers
#' @export
stan_binomial <- function(formula, 
                          data,
                          link = "logit",
                          ...,
                          prior = normal(),
                          prior_intercept = normal(),
                          prior_aux = cauchy(0, 5),
                          prior_covariance = decov(),
                          prior_PD = FALSE,
                          algorithm = c("sampling", "optimizing", 
                                        "meanfield", "fullrank"),
                          adapt_delta = NULL,
                          QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc$link <- NULL
  mc$family <- binomial(link = link)
  if (has_varying_terms(formula)) {
    mc[[1L]] <- quote(stan_glmer)
  } else {
    mc[[1L]] <- quote(stan_glm)
    mc$prior_covariance <- NULL
  }
  out <- eval(mc, parent.frame())
  out$call <- call
  return(out)
}

#' @rdname model-wrappers
#' @export
stan_gaussian <- function(formula, 
                              data,
                              link = "identity",
                              ...,
                              prior = normal(),
                              prior_intercept = normal(),
                              prior_aux = cauchy(0, 5),
                              prior_covariance = decov(),
                              prior_PD = FALSE,
                              algorithm = c("sampling", "optimizing", 
                                            "meanfield", "fullrank"),
                              adapt_delta = NULL,
                              QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc$link <- NULL
  mc$family <- gaussian(link = link)
  if (has_varying_terms(formula)) {
    mc[[1L]] <- quote(stan_glmer)
  } else {
    mc[[1L]] <- quote(stan_glm)
    mc$prior_covariance <- NULL
  }
  out <- eval(mc, parent.frame())
  out$call <- replace_function_in_call(call)
  return(out)
}



#' @rdname model-wrappers
#' @export
stan_gamma <- function(formula, 
                           data,
                             link = "inverse",
                             ...,
                             prior = normal(),
                             prior_intercept = normal(),
                             prior_aux = cauchy(0, 5),
                             prior_covariance = decov(),
                             prior_PD = FALSE,
                             algorithm = c("sampling", "optimizing", 
                                           "meanfield", "fullrank"),
                             adapt_delta = NULL,
                             QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc$link <- NULL
  mc$family <- Gamma(link = link)
  if (has_varying_terms(formula)) {
    mc[[1L]] <- quote(stan_glmer)
  } else {
    mc[[1L]] <- quote(stan_glm)
    mc$prior_covariance <- NULL
  }
  out <- eval(mc, parent.frame())
  out$call <- call
  return(out)
}

#' @rdname model-wrappers
#' @export
stan_inv_gaussian <- function(formula, 
                             data,
                             link = "1/mu^2",
                             ...,
                             prior = normal(),
                             prior_intercept = normal(),
                             prior_aux = cauchy(0, 5),
                             prior_covariance = decov(),
                             prior_PD = FALSE,
                             algorithm = c("sampling", "optimizing", 
                                           "meanfield", "fullrank"),
                             adapt_delta = NULL,
                             QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc$link <- NULL
  mc$family <- inverse.gaussian(link = link)
  if (has_varying_terms(formula)) {
    mc[[1L]] <- quote(stan_glmer)
  } else {
    mc[[1L]] <- quote(stan_glm)
    mc$prior_covariance <- NULL
  }
  out <- eval(mc, parent.frame())
  out$call <- call
  return(out)
}



#' @rdname model-wrappers
#' @export
stan_poisson <- function(formula, 
                             data,
                             link = "log",
                             ...,
                             prior = normal(),
                             prior_intercept = normal(),
                             prior_aux = cauchy(0, 5),
                             prior_covariance = decov(),
                             prior_PD = FALSE,
                             algorithm = c("sampling", "optimizing", 
                                           "meanfield", "fullrank"),
                             adapt_delta = NULL,
                             QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc$link <- NULL
  mc$family <- poisson(link = link)
  if (has_varying_terms(formula)) {
    mc[[1L]] <- quote(stan_glmer)
  } else {
    mc[[1L]] <- quote(stan_glm)
    mc$prior_covariance <- NULL
  }
  out <- eval(mc, parent.frame())
  out$call <- call
  return(out)
}

#' @rdname model-wrappers
#' @export
stan_neg_binomial <- function(formula, 
                        data,
                        link = "log",
                        ...,
                        prior = normal(),
                        prior_intercept = normal(),
                        prior_aux = cauchy(0, 5),
                        prior_covariance = decov(),
                        prior_PD = FALSE,
                        algorithm = c("sampling", "optimizing", 
                                      "meanfield", "fullrank"),
                        adapt_delta = NULL,
                        QR = FALSE) {
  if ("family" %in% names(list(...)))
    stop("'family' should not be specified.")
  mc <- call <- match.call()
  if (!"formula" %in% names(call))
    names(call)[2L] <- "formula"
  mc$link <- NULL
  mc$family <- neg_binomial_2(link = link)
  if (has_varying_terms(formula)) {
    mc[[1L]] <- quote(stan_glmer)
  } else {
    mc[[1L]] <- quote(stan_glm)
    mc$prior_covariance <- NULL
  }
  out <- eval(mc, parent.frame())
  out$call <- call
  return(out)
}



# internal ----------------------------------------------------------------
has_varying_terms <- function(formula) {
  bars <- lme4::findbars(formula)
  isTRUE(length(bars) > 0)
}

replace_function_in_call <- function(x) {
  x[[1]] <- if (has_varying_terms(x[["formula"]])) 
    quote(stan_glmer) else quote(stan_glm)
  x
}

