

fusedlasso <- function(x, y, lambda.lasso = 0, lambda.fused = 0, groups = NULL
                       family = c("gaussian", "binomial", "multinomial"), 
                       opts = NULL, class.weights = NULL) {
  family <- match.arg(family)
  colM <- colMeans(x)
  p <- ncol(x)
  n <- nrow(x)
  x.tilde <- x - matrix(rep(colM, n), ncol=p, byrow=TRUE)
  if (is.null(opts)) {
    opts <- sllOpts()
  }
  opts$fusedPenalty <- lambda.fused
  
  if (family == "gaussian") {
    res <- fusedLeastR(x = x.tilde, y = y, lambda = lambda.lasso, opts = opts)
    res$intercept <- mean(y)
  } else if (family == "binomial") {
    res <- fusedLogisticR(x = x.tilde, y, lambda = lambda.lasso, 
                          class.weights = class.weights, opts = opts)
  } else if (family == "multinomial") {
    res <- fusedMultinomialLogistic(x = x.tilde, y, lambda = lambda.lasso, 
                                    groups = groups, opts = opts)
  }
  res
}