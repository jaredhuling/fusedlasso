fusedMultinomialLogistic <- function(x, y, lambda, lambda.group = 0, groups = NULL, 
                                     class.weights = NULL, opts=NULL) {
  
  sz <- dim(x)
  n <- sz[1]
  p <- sz[2]
  
  y.f <- as.factor(y)
  #if (is.factor(y)) {
  #  y <- levels(y)[y]
  #}
  classes <- levels(y.f)
  K <- length(classes)
  
  y.mat <- array(-1, dim = c(n, K))
  for (k in 1:K) {
    y.mat[y.f == classes[k],k] <- 1
  }
  
  betas <- array(NA, dim = c(K, p))
  intercepts <- numeric(K)
  
  stopifnot(lambda > 0)
  #if ( any(sort(unique(y)) != c(-1, 1)) ) {
  #  stop("y must be in {-1, 1}")
  #}
  opts.orig <- opts
  
  
  # if groups are given, get unique groups
  
  if (!is.null(groups)) {
    unique.groups <- vector(mode = "list", length = K)
    if (is.list(groups)) {
      if (length(groups) != K) {
        stop("Group list but have one element per class")
      }
      
      for (k in 1:K) {
        unique.groups[[k]] <- sort(unique(groups[!is.na(groups[[k]])]))
      }
    } else {
      unique.groups[1:K] <- sort(unique(groups[!is.na(groups)]))
      gr.list <- vector(mode = "list", length = K)
      gr.list[1:K] <- rep(list(groups), K)
      groups <- gr.list
    }
  }

  # run sllOpts to set default values (flags)
  opts <- sllOpts(opts.orig)
  
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
  
  
  
  #code current y as 1 and -1
  #y.k <- 2 * (y.f == classes[k]) - 1
  
  ## Group & others
  
  # The parameter 'weight' contains the weight for each training sample.
  # See the definition of the problem above.
  # The summation of the weights for all the samples equals to 1.
  p.flag <- (y.mat == 1)
  if (!is.null(class.weights)) {
   sWeight <- class.weights
   
   if (length(class.weights) != 2 || class.weights[1] <= 0 || class.weights[2] <= 0) {
     stop("class weights must contain 2 positive values")
   }
   
   weight <- numeric(n)
   
   m1 <- sum(p.flag) * sWeight[1]
   m2 <- sum(!p.flag) * sWeight[2]
   weight[which(p.flag)] <- sWeight[1] / (m1 + m2)
   weight[which(!p.flag)] <- sWeight[2] / (m1 + m2)
  } else {
   weight <- array(1, dim = c(n, K)) / n
  }
  
  ## L2 norm regularization
  if (!is.null(opts$rsL2)) {
    rsL2 <- opts$rsL2
    if (rsL2 < 0) {
      stop("opts$rsL2 must be nonnegative")
    }
  } else {
    rsL2 <- 0
  }
  
  #m1 <- sum(weight[which(p.flag)])
  m1 <- colSums(p.flag * weight)
  m2 <- 1 - m1
  
  
  ## L1 norm regularization
  if (opts$rFlag != 0) {
    if (lambda < 0 || lambda > 1) {
      stop("opts.rFlag=1, and lambda should be in [0,1]")
    }
    
    ## we compute ATb for computing lambda_max, when the input lambda is a ratio
    b <- array(0, dim = c(n, K))
    b[which(p.flag)] <- m2
    b[which(!p.flag)] <- -m1
    b <- b * weight
    #b <- b / n
    
    ## compute xTb
    if (opts$nFlag == 0) {
      xTb <- crossprod(x, b)
    } else if (opts$nFlag == 1) {
      xTb <- (crossprod(x, b) - colSums(b) * mu) / nu
    } else {
      invNu <- b / nu
      xTb <- crossprod(x, invNu) - colSums(invNu) * mu
    }
    
    lambda.max <- max(abs(xTb))
    lambda <- lambda * lambda.max
    
    if (is.null(opts$fusedPenalty)) {
      lambda2 <- 0
    } else {
      lambda2 <- lambda.max * opts$fusedPenalty
    }
    
    rsL2 <- rsL2 * lambda.max
    
  } else {
    if (is.null(opts$fusedPenalty)) {
      lambda2 <- 0
    } else {
      lambda2 <- opts$fusedPenalty
    }
  }
  
  ## initialize a starting point
  if (opts$init == 2) {
    beta <- array(0, dim = c(p, K))
    c <- log(m1 / m2)
  } else {
    if (!is.null(opts$x0)) {
      beta <- opts$x0
      if (any(dim(beta) != c(p, K))) {
        stop("initialization must be of dimension p x K")
      }
    } else {
      beta <- array(0, dim = c(p, K))
    }
    
    if (!is.null(opts$c0)) {
      c <- opts$c0
      if (length(c) != K) {
        stop("c0 must be of length K")
      }
    } else {
      c <- log(m1 / m2)
    }
  }
  
  ## Compute x beta
  if (opts$nFlag == 0) {
    xbeta <- x %*% beta
  } else if (opts$nFlag == 1) {
    invNu <- beta / nu
    mu.invNu <- as.double(crossprod(mu, invNu))
    xbeta <- x %*% invNu - rep(mu.invNu, n)
  } else {
    mubeta <- as.double(crossprod(mu, beta))
    xbeta <- (x %*% beta - rep(mubeta, n)) / nu
  }
  
  # restart the program for better efficiency
  # this is a newly added function
  if (is.null(opts$rStartNum)) {
    opts$rStartNum <- opts$maxIter
  } else {
    if (opts$rStartNum <= 0) {
      opts$rStartNum <- opts$maxIter
    }
  }
  
  ### THE MAIN CODE
  
  # z0 is the starting point for flsa
  z0 <- array(0, dim = c((p-1), K))
  
  ValueL <- funVal <- numeric(opts$maxIter)
  
  ## The Armijo Goldstein line search scheme + accelerated gradient descent
  if (opts$mFlag == 0 && opts$lFlag == 0) {
    
    # this flag tests whether the gradient step only changes a little
    bFlag <- 0 
    
    # the intial guess of the Lipschitz continuous gradient
    L <- 1 / n + rsL2
    
    # the product between weight and y
    #weighty <- weight * y.k
    
    #initial assignments
    betap <- beta
    xbetap <- xbeta
    betabetap <- numeric(p)
    cp <- c; ccp <- numeric(K)
    
    alphap <- 0; alpha <- 1
    
    for (iterStep in 1:opts$maxIter) {
      
      diff.prev <- -n
      
      ## --------------------------- step 1 ---------------------------
      ##      compute search point s based on xp and x (with beta)
      bet <- (alphap - 1) / alpha
      s <- beta + bet * betabetap
      sc <- c + bet * ccp
      
      ## --------------------------- step 2 ---------------------------
      ##  line search for L and compute the new approximate solution x
      
      # compute xs = x * s
      xs <- xbeta + bet * (xbeta - xbetap)
      
      #aa <- -y.k * (xs + sc)
      aa <- (xs + rep(sc, each = n))
      
      # fun.s is the logistic loss at the search point
      bb <- pmax(- y.mat * aa, 0)
      #fun.s <- as.double( crossprod(weight, (log(exp(-bb) + exp(aa - bb)) + bb)) ) + 
      #  ( rsL2 / 2 ) * as.double(crossprod(s))

      
      # compute prob=[p_1;p_2;...;p_n]
      #prob <- 1 / (1 + exp(aa))
      prob <- exp(aa)
      

      rSp <- rowSums(prob)
      
      fun.s <- -sum(rowSums(((y.mat + 1) / 2) * aa) - log( rSp )) / n + 
        ( rsL2 / 2 ) * sum(as.double(crossprod(s)))
      
      prob <- prob / rSp

      #b <- -weighty * (1 - prob)
      b <- -((y.mat+1)/2 - prob) * weight

      #the gradient of c
      #gc <- (colSums(y.mat+1)/(2) - 1) / n
      gc <- colSums(b)

      #  should be sum i=1:n { sum k=1:K {y_i^(k)} - p_ij} 
      
      #compute g= xT b, the gradient of beta
      if (opts$nFlag == 0) {
        g <- crossprod(x, b)
      } else if (opts$nFlag == 1) {
        g <- (crossprod(x, b) - colSums(b) * mu) / nu
      } else {
        invNu <- b / nu
        g <- crossprod(x, invNu) - colSums(invNu) * mu
      }
      #add the squared L2 norm regularization term
      g <- g + rsL2 * s

      #assignments
      betap <- beta
      xbetap <- xbeta
      cp <- c
      
      while (TRUE) {
        # let s walk in a step in the antigradient of s to get v
        # and then do the Lq/L1-norm regularized projection
        v <- s - g / L
        c <- sc - gc / L

        for (k in 1:K) {
          if (is.null(groups)) {
            res <- flsa(v[, k], z0[, k], lambda / L, lambda2 / L, p,
                        1000, 1e-8, 1, 6)
            beta[, k] <- res[[1]]
            z0[, k] <- res[[2]]
            infor <- res[[3]]
          } else {
            
            if (any(is.na(groups[[k]]))) {
              ## don't apply fused lasso penalty
              ## to variables with group == NA 
              gr.idx <- which(is.na(groups[[k]]))
              gr.p <- length(gr.idx)
              if (any(gr.idx == 1)) {
                gr.idx.z <- gr.idx[gr.idx != 1] - 1
              } else {
                gr.idx.z <- gr.idx[-gr.p]
              }
              
              res <- flsa(v[gr.idx, k], z0[gr.idx.z, k], lambda / L, 0, gr.p,
                          1000, 1e-8, 1, 6)            
              
              beta[gr.idx, k] <- res[[1]]
              z0[gr.idx.z, k] <- res[[2]]
              infor <- res[[3]]
            }
            
            for (t in 1:length(unique.groups[[k]])) {
              gr.idx <- which(groups[[k]] == unique.groups[[k]][t])
              gr.p <- length(gr.idx)
              if (any(gr.idx == 1)) {
                gr.idx.z <- gr.idx[gr.idx != 1] - 1
              } else {
                gr.idx.z <- gr.idx[-gr.p]
              }
              
              res <- flsa(v[gr.idx, k], z0[gr.idx.z, k], lambda / L, lambda2 / L, gr.p,
                          1000, 1e-8, 1, 6)
              
              if (lambda.group > 0) {
                ## 2nd Projection:
                ## argmin_w { 0.5 \|w - w_1\|_2^2
                ##          + lambda_3 * \|w_1\|_2 }
                ## This is a simple thresholding:
                ##    w_2 = max(\|w_1\|_2 - \lambda_3, 0)/\|w_1\|_2 * w_1
                nm = norm(res[[1]], type = "2")
                if (nm == 0) {
                  newbeta = numeric(length(res[[1]]))
                } else {
                  #apply soft thresholding, adjust penalty for size of group
                  newbeta = pmax(nm - lambda.group * sqrt(gr.p), 0) / nm * res[[1]]
                }
                end
              } else {
                newbeta <- res[[1]]
              }
              
              beta[gr.idx, k] <- newbeta
              z0[gr.idx.z, k] <- res[[2]]
              infor <- res[[3]]
            }
          }
        } #end loop over classes

        # the difference between the new approximate 
        # solution x and the search point s
        v <- beta - s

        ## Compute x beta
        if (opts$nFlag == 0) {
          xbeta <- x %*% beta
        } else if (opts$nFlag == 1) {
          invNu <- beta / nu
          mu.invNu <- as.double(crossprod(mu, invNu))
          xbeta <- x %*% invNu - rep(mu.invNu, n)
        } else {
          mubeta <- as.double(crossprod(mu, beta))
          xbeta <- (x %*% beta - rep(mubeta, n)) / nu
        }
        
        #aa <- -y.k * (xbeta + c)
        aa <- (xbeta + rep(c, each = n))
        
        # fun.beta is the logistic loss at the new approximate solution
        bb <- pmax(- y.mat * aa, 0)

        
        fun.beta <- -sum(rowSums(((y.mat + 1) / 2) * aa) - log( rowSums(exp(aa)) )) / n + 
          ( rsL2 / 2 ) * sum(as.double(crossprod(beta)))
        

        
        r.sum <- norm(v, type = "F") ^ 2 / 2 + sum((c - sc)^2) / 2
        fzp.gamma <- fun.s + sum(sum(v * g)) + L * r.sum + sum((c - sc) * gc)


        if (r.sum <= 1e-18) {
          #this shows that the gradient step makes little improvement
          bFlag <- 1
          break
        }
        
        ## the condition is fun.beta <= fun.s + v'* g + c * gc
        ##                           + L/2 * (v'*v + (c-sc)^2 )
        
        if (fun.beta <= fzp.gamma ) { #| funval.s - funval < diff.prev
          break
        } else {
          #L <- max(2 * L, (fun.beta * L) / fzp.gamma)
          L <- 2 * L
        }
        
        #diff.prev <- funval.s - funval
        
      } # end while loop
      
      ## --------------------------- step 3 ---------------------------
      ##      update alpha and alphap, and check for convergence
      alphap <- alpha
      alpha <- (1 + sqrt(4 * alpha * alpha + 1)) / 2
      
      ## store values for L
      ValueL[iterStep] <- L
      
      betabetap <- beta - betap
      ccp <- c - cp
      
      fused.pen <- group.pen <- 0
      for (t in 1:length(unique.groups[[k]])) {
        gr.idx <- which(groups[[k]] == unique.groups[[k]][t])
        gr.p <- length(gr.idx)
        if (gr.p > 1) {
          fused.pen <- fused.pen + sum(abs(beta[gr.idx[2:(gr.p)], k] - beta[gr.idx[1:(gr.p - 1)], k]))
          group.pen <- group.pen + sqrt(sum(beta[gr.idx, k] ^ 2))
        }
      }
      
      funVal[iterStep] <- fun.beta + lambda * sum(abs(beta)) +
        lambda2 * fused.pen + lambda.group * group.pen
      
      if (bFlag) {
        break
      }
      
      tf <- opts$tFlag
      
      if (tf == 0) {
        if (iterStep > 2) {
          if (abs( funVal[iterStep] - funVal[iterStep - 1] ) <= opts$tol) {
            break
          }
        }
      } else if (tf == 1) {
        if (iterStep > 2) {
          if (abs( funVal[iterStep] - funVal[iterStep - 1] ) <= 
                opts$tol * funVal[iterStep - 1]) {
            break
          }
        }
      } else if (tf == 2) {
        if (funVal[iterStep] <= opts$tol) {
          break
        }
      } else if (tf == 3) {
        norm.bbp <- sqrt(as.double(crossprod(bbp)))
        if (norm.bbp <= opts$tol) {
          break
        }
      } else if (tf == 4) {
        norm.bp <- sqrt(as.double(crossprod(bp)))
        norm.bbp <- sqrt(as.double(crossprod(bbp)))
        if (norm.bbp <= opts$tol * max(norm.bp)) {
          break
        }
      } else if (tf == 5) {
        if (iterStep >= opts$maxIter) {
          break
        }
      }
      
      # restart the program every opts$rStartNum
      if ((iterStep %% opts$rStartNum == 0)) {
        alphap <- 0; alpha <- 1
        betap <- beta; xbetap <- xbeta; L <- 1
        betabetap <- numeric(p)
      }
      
    } #end of iterStep loop
    
    funVal <- funVal[1:iterStep]
    ValueL <- ValueL[1:iterStep]
    
  } else {
    stop("This function only supports opts.mFlag=0 & opts.lFlag=0")
  }

  
  list(beta = t(beta), intercept = c, funVal = funVal, ValueL = ValueL, bFlag = bFlag)
}