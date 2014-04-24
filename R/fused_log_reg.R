
fusedLassoLogistic <- function(x, y, lambda, class.weights = NULL, opts=NULL) {
  
  sz <- dim(x)
  n <- sz[1]
  p <- sz[2]
  
  stopifnot(lambda > 0)
  if ( any(sort(unique(y)) != c(-1, 1)) ) {
    stop("y must be in {-1, 1}")
  }
  
  # run sllOpts to set default values (flags)
  opts <- sllOpts(opts)
  
  ## Set up options
  if (opts$nFlag != 0) {
    if (!is.null(opts$mu)) {
      mu <- opts$mu
      stopifnot(length(mu) == p)
    } else {
      mu <- colMeans(x)
    }
    
    if (opts$nFlag == 1) {
      if (!is.null(opts$nu)) {
        nu <- opts$nu
        stopifnot(length(nu) == p)
      } else {
        nu <- sqrt(colSums(x^2) / n)
      }
    }
    
    if (opts$nFlag == 2) {
      if (!is.null(opts$nu)) {
        nu <- opts$nu
        stopifnot(length(nu) == n)
      } else {
        nu <- sqrt(rowSums(x^2) / p)
      }
    }
    
    ## If some values of nu are small, it might
    ## be that the entries in a given row or col
    ## of x are all close to zero. For numerical 
    ## stability, we set these values to 1.
    ind.zero <- which(abs(nu) <= 1e-10)
    nu[ind.zero] <- 1
    
  }
  
  ## Group & others
  
  # The parameter 'weight' contains the weight for each training sample.
  # See the definition of the problem above.
  # The summation of the weights for all the samples equals to 1.
  if (!is.null(class.weights)) {
    sWeight <- class.weights
    
    if (length(class.weights) != 2 || class.weights[1] <= 0 || class.weights[2] <= 0) {
      stop("class weights must contain 2 positive values")
    }
    
    weight <- numeric(n)
    p.flag <- (y == 1)
    m1 <- sum(p.flag) * sWeight[1]
    m2 <- sum(!p.flag) * sWeight[2]
    weight[which(p.flag)] <- sWeight[1] / (m1 + m2)
    weight[which(!p.flag)] <- sWeight[2] / (m1 + m2)
  } else {
    weight <- rep(1, n) / n
  }
  
  
  
}