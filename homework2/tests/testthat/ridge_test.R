library(testthat)
library(MASS)

context("Test the output of ridge_regression().")

test_that("Your ridge_regression() function works in an easy case.", {

  data(iris)

  fit_lm_ridge <- lm.ridge(Sepal.Length ~ ., iris, lambda = 0.1)

  fit_ridge_regression <- ridge_regression(Sepal.Length  ~ ., iris, lambda = 0.1)

  expect_equivalent(coef(fit_lm_ridge), fit_ridge_regression,
                    tolerance = 1e-3)
})

test_that("Your ridge_regression() function works with contrasts.", {

  data(iris)

  fit_lm_ridge <- MASS::lm.ridge(Sepal.Length ~ ., iris,
                                 contrasts = list(Species = "contr.sum"))

  fit_ridge_regression <- ridge_regression(Sepal.Length  ~ ., iris, contrasts = list(Species = "contr.sum"))

  expect_equivalent(coef(fit_lm_ridge), fit_ridge_regression,
                    tolerance = 1e-3)

})
