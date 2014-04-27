
#' Efficient Fused Lasso Algorithm (EFLA) 
#'
#' @param x input matrix. Each row is an observation, each column corresponds to a covariate
#' @param y response vector of length nobs. numeric if \code{family == "gaussian"}, vector in \{-1, 1\} if \code{family == "binomial"}, 
#' vector with values in \{1, 2, ..., k\} if \code{family == "multinomial"}
#' @param lambda.lasso tuning parameter for lasso penalty
#' @param lambda.fused tuning parameter for fused lasso penalty
#' @param groups vector which defines the grouping of the variables for which to apply the fused lasso penalty. Components sharing the
#' same number build a group. Non-fused-lasso-penalized coefficients are marked with NA. 
#' Currently only works for \code{family == "multinomial"}
#' @param family "gaussian" for linear regression, "binomial" for logistic regression, "multinomial" for multinomial logistic regression
#' @param opts options as defined by \code{sllOpts()} function
#' @return 
#' \itemize{
#'  \item{\code{beta}} {p vector (or K x p matrix) of estimated coefficients}
#'  \item{\code{intercept}} {estimated intercept}
#' }
#' @references \emph{An Efficient Algorithm for a Class of Fused Lasso Problems}, Liu et al. 2010 
#' http://www.public.asu.edu/~jye02/Publications/Papers/rp589f-liu.pdf
#' 
#' C code by Jun Liu, Shuiwang Ji, and Jieping Ye
#' http://www.public.asu.edu/~jye02/Software/SLEP/index.htm
#' @export
#' @examples
#' nobs <- 10000
#' nvars <- 50
#' 
#' #generate data
#' set.seed(123)
#' true.beta <- rnorm(nvars) * rbinom(nvars, 1, 0.25)
#' true.beta[1:3] <- c(1, 1.06, 0.95)
#' x <- matrix(rnorm(nobs * nvars), ncol = nvars)
#' 
#' #generate binary outcome
#' log.p.ratio <- x %*% true.beta
#' prob.y.1 <- 1 / (1 + exp(-log.p.ratio))
#' y <- rbinom(nobs, 1, prob = prob.y.1)
#' y1 <- ifelse(y == 0, -1, y)
#' 
#' #fit fused lasso logistic model
#' res <- fusedlasso(x, y1, lambda.lasso = 0.005, lambda.fused = 0.01, family = "binomial")
#' round(res$beta, 5)
fusedlasso <- function(x, y, weights = rep(1, nrow(x)), 
                       lambda.lasso = 0, lambda.fused = 0, groups = NULL,
                       family = c("gaussian", "binomial", "multinomial"), 
                       opts = NULL, class.weights = NULL) {
  
  family <- match.arg(family)
  colM <- colMeans(sqrt(weights) * x)
  p <- ncol(x)
  n <- nrow(x)

  x.tilde <- sqrt(weights) * x - matrix(rep(colM, n), ncol=p, byrow=TRUE)
  if (is.null(opts)) {
    opts <- sllOpts()
  }
  opts$fusedPenalty <- lambda.fused
  
  if (family == "gaussian") {
    res <- fusedLeastR(x = x.tilde, y = sqrt(weights) * y, 
                       lambda = lambda.lasso, 
                       groups = groups, opts = opts)
    res$intercept <- mean(sqrt(weights) * y)
  } else if (family == "binomial") {
    
    if (any(sort(unique(y)) != c(-1, 1) )) {
      stop ("y must take values {-1, 1}")
    }
    
    res <- fusedLogisticR(x = x.tilde, y, lambda = lambda.lasso, 
                          class.weights = class.weights, opts = opts)
  } else if (family == "multinomial") {
    res <- fusedMultinomialLogistic(x = x.tilde, y, lambda = lambda.lasso, 
                                    groups = groups, opts = opts)
  }
  res
}