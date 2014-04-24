

initFactor <- function(b.norm, xb, y, z, 
                       funName = c("LeastC", "LeastR", "glLeastR",
                                   "mcLeastR", "mtLeastR", "nnLeastR",
                                   "nnLeastC", "mcLeastC"),
                       rsL2, b.2norm) {
  
  ## function initFactor
  #     compute the an optimal constant factor for the initialization
  #
  #
  # Input parameters:
  # x_norm-      the norm of the starting point
  # Ax-          A*x, with x being the initialization point
  # y-           the response matrix
  # z-           the regularization parameter or the ball
  # funName-     the name of the function
  #
  # Output parameter:
  # ratio-       the computed optimal initialization point is ratio*x
  #
  ## Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
  #
  # For any problem, please contact with Jun Liu via j.liu@asu.edu
  #
  # Last revised on August 2, 2009.
  
  funName <- match.arg(funName)
  
  if (funName == "LeastC") {
    ratio.max <- z / b.norm
    ratio.optimal <- as.double(crossprod(xb, y)) / 
      (as.double(crossprod(xb)) + rsL2 * b.2norm )
    
    if (abs(ratio.optimal) <= ratio.max) {
      ratio <- ratio.optimal
    } else if(ratio.optimal < 0) {
      ratio <- -ratio.max
    } else {
      ratio <- ratio.max
    }
    
  } else if (funName == "LeastR") {
    ratio <- (as.double(crossprod(xb, y)) - z * b.norm) /
      (as.double(crossprod(xb)) + rsL2 * b.2norm )
  } else if (funName == "glLeastR") {
    ratio <- (as.double(crossprod(xb, y)) - z * b.norm) /
      (as.double(crossprod(xb)) )
  } else {
    cat(funName, "not supported yet")
    stop(".")
  }
  ratio
}