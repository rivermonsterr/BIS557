#' Ridge Regression
#'
#'@description an implementation of ridge regression using singular value decomposition
#' @param form formula y ~ .
#' @param data dataframe
#' @param lambda ridge penalty term; default = 0
#' @param contrasts a list of contrasts
#'
#' @return List of beta coefficients
#' @export
#'

ridge_regression <- function(formula, data, lambda = 0, contrasts = NULL) {
  # mitigate issue with folds
  rownames(data) <- NULL

  # model matrix
  X <- model.matrix(formula, data, contrasts.arg= contrasts)

  # response variable
  Y <- data[[as.character(formula)[2]]][as.numeric(rownames(X))]

  # rescaling
  Xmean <- colMeans(X[, -1])
  Ymean <- mean(Y)

  # match lm.ridge function
  X <- X[, -1] - rep(Xmean, rep(nrow(X), ncol(X) - 1))
  Xscale <- drop(rep(1/nrow(X), nrow(X)) %*% X^2)^0.5
  X <- X/rep(Xscale, rep(nrow(X), ncol(X)))

  Y <- Y - Ymean

  # beta matrix

  coef <- matrix(NA_real_, nrow = length(lambda), ncol= ncol(X))

  # svd to determine coefficients

  svd<- svd(X)
  D<- diag(svd$d / (svd$d^2 + lambda))
  coef<- svd$v %*% D %*% t(svd$u) %*% Y

  scaledcoef <- t(as.matrix(coef / Xscale))

  # add intercept
  intercept <- Ymean - scaledcoef %*% Xmean
  coef <- cbind(intercept, scaledcoef)

  # create output
  coef<- as.vector(coef)
  names(coef) <- c("Intercept", colnames(X))
  attributes(coef)$formula <- formula
  class(coef)<- c(class(coef), "ridge_regression")
  coef

}

predict.ridge_regression <- function(object, ...) {
  dots <- list(...)
  x_frame <- dots[[1]]
  if (!is.data.frame(x_frame)) {
    stop(red("The first argument should be a data.frame of values",
             "to predict"))
  }
  X <- model.matrix(attributes(object)$formula, x_frame)
  X %*% object
}
