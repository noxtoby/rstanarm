# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2017 Sam Brilleman
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
library(lme4)
ITER <- 1000
CHAINS <- 1
SEED <- 12345
REFRESH <- 0L
set.seed(SEED)
if (interactive()) 
  options(mc.cores = parallel::detectCores())

TOLSCALES <- list(
  lmer_fixef = 0.25,  # how many SEs can stan_jm fixefs be from lmer fixefs
  lmer_ranef = 0.05, # how many SDs can stan_jm ranefs be from lmer ranefs
  glmer_fixef = 0.3, # how many SEs can stan_jm fixefs be from glmer fixefs
  glmer_ranef = 0.1 # how many SDs can stan_jm ranefs be from glmer ranefs
)

source(test_path("helpers", "expect_matrix.R"))
source(test_path("helpers", "expect_stanreg.R"))
source(test_path("helpers", "expect_stanmvreg.R"))
source(test_path("helpers", "expect_survfit.R"))
source(test_path("helpers", "expect_ppd.R"))
source(test_path("helpers", "expect_identical_sorted_stanmats.R"))
source(test_path("helpers", "SW.R"))
source(test_path("helpers", "get_tols.R"))
source(test_path("helpers", "recover_pars.R"))

context("stan_ltjmm")

#----  Data (for non-Gaussian families)

pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
pbcLong$ybino <- as.integer(rpois(nrow(pbcLong), 5))
pbcLong$ypois <- as.integer(pbcLong$albumin)
pbcLong$ynbin <- as.integer(rnbinom(nrow(pbcLong), 3, .3))
pbcLong$ygamm <- as.numeric(pbcLong$platelet / 10)
pbcLong$xbern <- as.numeric(pbcLong$platelet / 100)
pbcLong$xpois <- as.numeric(pbcLong$platelet / 100)
pbcLong$xgamm <- as.numeric(pbcLong$logBili)

#----  Models

# univariate LTJMM
fm1 <- logBili ~ I(year+lt) + (year | id)
o<-SW(m1 <- stan_ltjmm(fm1, pbcLong, lt_var='lt', lt_formula = ~I(year+lt), 
  iter = 10, chains = 1, seed = SEED))

# multivariate LTJMM
fm2 <- list(logBili ~  I(year+lt) + (year | id), albumin ~  I(year+lt) + (year | id))
o<-SW(m2 <- stan_ltjmm(fm2, pbcLong, lt_var='lt', lt_formula = ~I(year+lt), 
  iter = 10, chains = 1, seed = SEED))

#----  Tests for stan_mvmer arguments

test_that("formula argument works", {
  SW(m991 <- update(m1, formula. = list(fm1)))
  expect_identical(as.matrix(m1), as.matrix(m991)) # fm as list
})

test_that("data argument works", {
  SW(m991 <- update(m1, data = list(pbcLong)))
  SW(m992 <- update(m2, data = list(pbcLong, pbcLong)))
  expect_identical(as.matrix(m1), as.matrix(m991)) # data as list
  expect_identical(as.matrix(m2), as.matrix(m992))
})

test_that("family argument works", {
  
  expect_output(ret <- update(m1, family = "gaussian"))
  expect_output(ret <- update(m1, family = gaussian))
  expect_output(ret <- update(m1, family = gaussian(link = identity)))
  
  expect_output(ret <- update(m1, formula. = ybern ~ ., family = binomial))
  expect_output(ret <- update(m1, formula. = ypois ~ ., family = poisson))
  expect_output(ret <- update(m1, formula. = ypois ~ ., family = neg_binomial_2))
  expect_output(ret <- update(m1, formula. = ygamm ~ ., family = Gamma, init = 0))
  expect_output(ret <- update(m1, formula. = ygamm ~ ., family = inverse.gaussian, init = 0))
  
  expect_error(ret <- update(m1, formula. = ybino ~ ., family = binomial))
  
  # multivariate model with combinations of family
  expect_output(ret <- update(m2, formula. = list(~ ., ybern ~ .), 
                              family = list(gaussian, binomial)))
})

test_that("prior_PD argument works", {
  expect_output(update(m1, prior_PD = TRUE))
})

test_that("adapt_delta argument works", {
  expect_output(update(m1, adapt_delta = NULL))
  expect_output(update(m1, adapt_delta = 0.8))
})

test_that("error message occurs for arguments not implemented", {
  expect_error(update(m1, weights = 1:10), "not yet implemented")
  expect_error(update(m1, QR = TRUE), "not yet implemented")
  expect_error(update(m1, sparse = TRUE), "not yet implemented")
})

#----  Compare estimates: univariate stan_ltjmm vs stan_glmer

if (interactive()) {
  test_that("coefs same for stan_glmer and stan_ltjmm", {
    SW(y1 <- stan_glmer(logBili ~ scale(year) + (1 | id), pbcLong, iter = 1000, 
      chains = CHAINS, seed = SEED))
    SW(y2 <- stan_ltjmm(logBili ~ I(scale(year)+lt) + (1 | id), pbcLong, iter = 1000, 
      chains = CHAINS, seed = SEED, lt_var='lt', lt_formula = ~I(scale(year)+lt)))
    expect_equivalent(fixef(y1), fixef(y2)$y1, tol=0.4)
    expect_equal(colMeans(log_lik(y1)), colMeans(log_lik(y2)), tol = 5)
  })
  test_that("coefs same for stan_glmer and stan_ltjmm, poisson", {
    SW(y1 <- stan_glmer(ypois ~ scale(year) + xpois + (1 | id), pbcLong, poisson, 
      init = 0, iter = 1000, chains = CHAINS, seed = SEED))
    SW(y2 <- stan_ltjmm(ypois ~ I(scale(year)+lt) + xpois + (1 | id), pbcLong, family=poisson, 
      iter = 1000, chains = CHAINS, seed = SEED, lt_var='lt', lt_formula = ~I(scale(year)+lt)))
    expect_equivalent(fixef(y1), fixef(y2)$y1, tol=0.03)
    expect_equal(colMeans(log_lik(y1)), colMeans(log_lik(y2)), tol = 0.03)
  })
  test_that("coefs same for stan_glmer and stan_ltjmm, negative binomial", {
    SW(y1 <- stan_glmer(ynbin ~ scale(year) + xpois + (1 | id), pbcLong, neg_binomial_2, 
      init = 0, iter = 1000, chains = CHAINS, seed = SEED))
    SW(y2 <- stan_ltjmm(ynbin ~ I(scale(year)+lt) + xpois + (1 | id), pbcLong, family=neg_binomial_2, 
      iter = 1000, chains = CHAINS, seed = SEED, lt_var='lt', lt_formula = ~I(scale(year)+lt)))
    expect_equivalent(fixef(y1), fixef(y2)$y1, tol=0.04)
    expect_equal(colMeans(log_lik(y1)), colMeans(log_lik(y2)), tol = 0.02)
  })
  test_that("coefs same for stan_jm and stan_glmer, Gamma", {
    # compare_glmer(ygamm ~ year + xgamm + (1 | id), Gamma(log))
    SW(y1 <- stan_glmer(ygamm ~ scale(year) + xgamm + (1 | id), pbcLong, Gamma(log), 
      init = 0, iter = 1000, chains = CHAINS, seed = SEED))
    SW(y2 <- stan_ltjmm(ygamm ~ I(scale(year)+lt) + xgamm + (1 | id), pbcLong, family=Gamma(log), 
      iter = 1000, chains = CHAINS, seed = SEED, lt_var='lt', lt_formula = ~I(scale(year)+lt)))
    expect_equivalent(fixef(y1), fixef(y2)$y1, tol=0.1)
    expect_equal(colMeans(log_lik(y1)), colMeans(log_lik(y2)), tol = 0.5)
  })
}

#----  Check methods and post-estimation functions

tmpdat <- pbcLong

o<-SW(f1 <- update(m1, formula. = list(logBili ~ I(year+lt) + (year | id)), data = tmpdat))

for (j in 1) {
  mod <- get(paste0("f", j))
  cat("Checking model:", paste0("f", j), "\n")

  expect_error(posterior_traj(mod), "stanjm")
  expect_error(posterior_survfit(mod), "stanjm")
     
  test_that("posterior_predict works with estimation data", {
    pp <- posterior_predict(mod, m = 1)
    expect_ppd(pp)
    if (mod$n_markers > 1L) {
      pp <- posterior_predict(mod, m = 2)
      expect_ppd(pp)
    }
  })  
  test_that("log_lik works with estimation data", {
    ll <- log_lik(mod)
    expect_matrix(ll)
    expect_identical(ll, log_lik(mod, m = 1))
    if (mod$n_markers > 1L)
      expect_matrix(log_lik(mod, m = 2))
  })   
  
  nd <- tmpdat[tmpdat$id == 2,]
  nd$lt <- 0
  test_that("posterior_predict works with new data (one individual)", {
    pp <- posterior_predict(mod, m = 1, newdata = nd)
    expect_ppd(pp)
    if (mod$n_markers > 1L) {
      pp <- posterior_predict(mod, m = 2, newdata = nd)
      expect_ppd(pp)
    }
  })     
  test_that("log_lik works with new data (one individual)", {
    ll <- log_lik(mod, newdata = nd)
    expect_matrix(ll)
    expect_identical(ll, log_lik(mod, m = 1, newdata = nd))
    if (mod$n_markers > 1L)
      expect_matrix(log_lik(mod, m = 2, newdata = nd))
    # log_lik is only designed for one submodel at a time so passing
    # newdata as a list should generate an error in validate_newdata
    expect_error(log_lik(mod, newdata = list(nd)), "data frame") 
  }) 
  
  nd <- tmpdat[tmpdat$id %in% c(1,2),]
  nd$lt <- 0
  test_that("posterior_predict works with new data (multiple individuals)", {
    pp <- posterior_predict(mod, m = 1, newdata = nd)
    expect_ppd(pp)
    if (mod$n_markers > 1L) {
      pp <- posterior_predict(mod, m = 2, newdata = nd)
      expect_ppd(pp)
    }
  })
  test_that("log_lik works with estimation data", {
    expect_matrix(log_lik(mod, newdata = nd))
    if (mod$n_markers > 1L)
      expect_matrix(log_lik(mod, m = 2, newdata = nd))
  }) 
  
  test_that("loo and waic work", {
    l <- suppressWarnings(loo(mod))
    w <- suppressWarnings(waic(mod))
    expect_s3_class(l, "loo")
    expect_s3_class(w, "loo")
    expect_s3_class(w, "waic")
    att_names <- c('names', 'dims', 'class', 'name', 'discrete', 'yhash', 'formula')
    expect_named(attributes(l), att_names)
    expect_named(attributes(w), att_names)
  })
  
  test_that("extraction methods work", {
    M <- mod$n_markers
    fe <- fixef(mod)
    re <- ranef(mod)
    ce <- coef(mod)
    mf <- model.frame(mod)
    tt <- terms(mod)
    fm <- formula(mod)
    fam <- family(mod)
    sig <- sigma(mod)
    expect_is(fe, "list"); expect_identical(length(fe), M)
    expect_is(re, "list"); expect_identical(length(re), M)
    expect_is(ce, "list"); expect_identical(length(re), M)
    expect_is(mf, "list"); expect_identical(length(mf), M); lapply(mf, function(x) expect_is(x, "data.frame"))
    expect_is(tt, "list"); expect_identical(length(tt), M); lapply(tt, function(x) expect_is(x, "terms"))
    expect_is(fm, "list"); expect_identical(length(fm), M); lapply(fm, function(x) expect_is(x, "formula"))
    expect_is(fam,"list"); expect_identical(length(fam),M); lapply(fam, function(x) expect_is(x, "family"))
    expect_is(sig, "numeric");
  })
  
  test_that("these extraction methods are currently disallowed", {
    expect_error(se(mod), "Not currently implemented")
    expect_error(fitted(mod), "Not currently implemented")
    expect_error(residuals(mod), "Not currently implemented")
  })
}

