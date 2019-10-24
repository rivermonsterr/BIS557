#' Cross-Validation for Ridge Regression
#'
#'@description methodology to find optimal lambdas using a given data set with ridge regression
#'
#' @param formula form "y~..."
#' @param data a data frame
#' @param folds the number of folds used for cross-validation
#' @param lambdas a vector of lambdas
#' @param contrasts
#'
#' @return a tibble containing a summary of the statistics
#' @importFrom stats sd model.matrix predict var
#' @importFrom rsample vfold_cv testing training
#' @importFrom foreach %dopar% %do%

#' @importFrom magrittr %>%
#' @import casl dplyr ggplot2
#' @export
#'
#' @examples
cross_validation <- function(formula, data, folds= 5, lambdas= seq(0, 0.5, 0.01), contrasts= NULL){
  # R CMD check
  i <- lambda <- `.` <- lower <- upper <- NULL

  # create folds
  folds <- vfold_cv(data,folds)

  # nested loop to find root mean squared error for each lambda for each fold

  rmses<- foreach(lambda = lambdas, .combine = rbind) %dopar% {
    foreach(i = seq_len(nrow(folds)), .combine = c) %do% {
      casl_util_rmse(
        testing(folds$splits[[i]])[[as.character(formula[2])]],
        predict(ridge_regression(formula, training(folds$splits[[i]]),
                                 lambda = lambda, contrasts= contrasts),
                testing(folds$splits[[i]])))
    }
  }

  # tibble of results
  results <- tibble(mean =apply(rmse, 1, mean),
                    sd = apply(rmse, 1, sd),
                    lambda = lambdas) %>%
    mutate(upper = mean + 2 * sd / nrow(.),
           lower = mean - 2 * sd / nrow(.))

  # lambda min
  lambda_min <- results$lambda[which.min(results$mean)]


  # find closest lambda 1se to the right
  find1se <- which((edf$mean > min(edf$mean) + sqrt(var(edf$mean))/length(edf$mean)))
  nextindex <- find1se[find1se > which.min(edf$mean)]
  lambda_1se <- edf$lambda[ifelse(!is.null(nextindex),
                                  min(find1se[find1se > which.min(edf$mean)]),
                                  which.min(edf$mean))]

  # create a list of output with table, lambda and rmse plot, lambda_min, and lambda_1se
  list(cv_table = edf, cv_plot = {
    ggplot2::ggplot(edf, aes(x = lambdas, y = mean, ymin = lower, ymax = upper)) +
      theme_minimal() +
      geom_errorbar(alpha = 0.3, width = 0) +
      geom_point(aes(color = "red")) +
      geom_vline(xintercept = lambda_min, linetype="dotted") +
      geom_vline(xintercept = lambda_1se, linetype="dotted") +
      ylab("Mean Squared Error") +
      xlab(expression(lambda)) +
      theme(legend.position = "none")
  },
  lambda_min = lambda_min,
  lambda_1se = lambda_1se)
}

