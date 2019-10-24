library(testthat)
library(glmnet)
library(Metrics)

context("Test the output of cross_validation().")

test_that("Your cross_validation() function finds the same lambda.min with the same inputs in an easy case.", {

  data(iris)

  fit_cv_glmnet <- cv.glmnet(model.matrix(Sepal.Length ~ ., iris), as.matrix(iris[,1]), alpha = 0, lambda = seq(0, 1, 0.03))

  fit_cv_ridge_regression <- cross_validation(Sepal.Length  ~ ., iris, folds = 2)

  expect_equivalent(fit_cv_glmnet$lambda.min, fit_cv_ridge_regression$lambda_min,
                    tolerance = 1e-3)
})

test_that("Your cross_validation() function works with contrasts.", {

  data(iris)

  fit_cv_glmnet <- cv.glmnet(model.matrix(Sepal.Length ~ ., iris,
                                          contrasts.arg = list(Species = "contr.sum")),
                             as.matrix(iris[,1]), alpha = 0, lambda = seq(0, 1, 0.03))

  fit_cv_ridge_regression <- cross_validation(Sepal.Length  ~ ., iris, contrasts = list(Species = "contr.sum"), folds = 2)

  expect_equivalent(fit_cv_glmnet$lambda.min, fit_cv_ridge_regression$lambda_min,
                    tolerance = 1e-3)
})
