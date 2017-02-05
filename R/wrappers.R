#' Create lightweight model wrappers
#' 
#' @name model-wrappers
#' 
#' @param family A \link[stats]{family} object, e.g.,
#'   \code{binomial(link="logit")}.
#' 
#' @return A modified version of one of the \pkg{rstanarm} modeling functions
#' (e.g. \code{stan_glm}) that uses the specified family.
#' 
#' @seealso \code{\link{stan_glm}}, \code{\link{stan_glmer}}, 
#' \code{\link{stan_gamm4}}
#' 
#' @examples
#' head(wells)
#' wells$dist100 <- wells$dist / 100
#' 
#' stan_logit <- stan_glm_wrapper(family = binomial(link = "logit"))
#' fit <- stan_logit(switch ~ dist100 + arsenic, data = wells, 
#'                   chains = 2, iter = 300) # for speed of example only
#' print(fit)
#' 
NULL

#' @rdname model-wrappers
#' @export
stan_glm_wrapper <- function(family) {
  .wrapper(rstanarm::stan_glm, family)
}

#' @rdname model-wrappers
#' @export
stan_glmer_wrapper <- function(family) {
  .wrapper(rstanarm::stan_glmer, family)
}

#' @rdname model-wrappers
#' @export
stan_gamm4_wrapper <- function(family) {
  .wrapper(rstanarm::stan_gamm4, family)
}


# internal ----------------------------------------------------------------
.wrapper <- function(fun, family) {
  stopifnot(is(family, "family"))
  wrapper <- match.fun(fun)
  formals(wrapper)[["family"]] <- family
  return(wrapper)
}
