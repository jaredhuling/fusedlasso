

fusedMultinomLogReg <- function(x, y, groups = NULL,
                                lambda.lasso = 0, lambda.fused = 0,
                                intercept = TRUE,
                                irls.maxiter = 30, irls.tol = 1e-10, 
                                beta.init = NULL, groups.in = NULL) {
  
  y.f <- as.factor(y)
  classes <- levels(y.f)
  K <- length(classes)
  
  nobs <- nrow(x)
  nvars <- ncol(x)
  len <- if (intercept) {nvars + 1} else {nvars}
  if (!is.null(beta.init)) {
    stopifnot(all(dim(beta.init) == c(K, len)))
  }
  w <- rep(0.5, nobs)
  betas <- if(is.null(beta.init)) {array(1, dim = c(K, len))} else {beta.init}
  beta <- betas[1,]
  converged <- rep(FALSE, K)
  for (i in 1:irls.maxiter) {
    prev <- betas
    for (k in 1:K) {
      if (!converged[k]) { 
        y.working <- 1 * (y.f == classes[k])
        
        init <- if (intercept) {prev[k,-1]} else {prev[k,]}
        
        beta.tmp <- fusedlasso(x, y.working, w, groups = groups,
                               lambda.lasso = lambda.lasso, 
                               lambda.fused = lambda.fused, 
                               family = "gaussian")
        
        if (intercept) {
          beta[-1] <- beta.tmp$beta
        } else {
          beta <- beta.tmp$beta
        }
        
        if (intercept) {
          xwb.tmp <- drop(x %*% beta[-1])
          #beta[1] <- mean( y.working - xwb.tmp)
          beta[1] <- beta.tmp$intercept
          xwb <- xwb.tmp + beta[1]
        } else {
          xwb <- drop(x %*% beta)
        }
        
        # update weights
        p <- 1 / (1 + exp(-xwb))
        w <- p * (1 - p)
        
        y.working <- xwb + (y - p) / w
        
        betas[k,] <- beta
        
        if (all(abs(beta - prev[k,]) < fused.tol)) {
          converged[k] <- TRUE
        }
        
      }
    }
    if (all(converged)) {
      cat("IRLS Converged at iteration: ", i, "\n")
      break
    }
    cat("IRLS iter: ", i, "\n")
  }
  betas
}