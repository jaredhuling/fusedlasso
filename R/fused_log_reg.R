
fusedLogisticR <- function(x, y, lambda, 
                           lambda.group = 0, groups = NULL,
                           class.weights = NULL, opts=NULL) {
  
  sz <- dim(x)
  n <- sz[1]
  p <- sz[2]
  
  stopifnot(lambda > 0)
  if ( any(sort(unique(y)) != c(-1, 1)) ) {
    stop("y must be in {-1, 1}")
  }
  
  
  
  # if groups are given, get unique groups
  if (!is.null(groups)) {
    unique.groups <- sort(unique(groups[!is.na(groups)]))
  }
  
  # run sllOpts to set default values (flags)
  opts <- sllOpts(opts)
  if (lambda.group > 0) {
    opts$tol <- 1e-10
  }
  
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
  p.flag <- (y == 1)
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
    weight <- rep(1, n) / n
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
  
  m1 <- sum(weight[which(p.flag)])
  m2 <- 1 - m1
  
  ## L1 norm regularization
  if (opts$rFlag != 0) {
    if (lambda < 0 || lambda > 1) {
      stop("opts.rFlag=1, and lambda should be in [0,1]")
    }
    
    ## we compute ATb for computing lambda_max, when the input lambda is a ratio
    b <- numeric(n)
    b[which(p.flag)] <- m2
    b[which(!p.flag)] <- -m1
    b <- b * weight
    
    ## compute xTb
    if (opts$nFlag == 0) {
      xTb <- crossprod(x, b)
    } else if (opts$nFlag == 1) {
      xTb <- (crossprod(x, b) - sum(b) * mu) / nu
    } else {
      invNu <- b / nu
      xTb <- crossprod(x, invNu) - sum(invNu) * mu
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
    beta <- numeric(p)
    c <- log(m1 / m2)
  } else {
    if (!is.null(opts$x0)) {
      beta <- opts$x0
      if (length(beta) != p) {
        stop("initialization must be of length p")
      }
    } else {
      beta <- numeric(p)
    }
    
    if (!is.null(opts$c0)) {
      c <- opts$c0
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
  z0 <- numeric(p-1)
  
  ValueL <- funVal <- numeric(opts$maxIter)
  
  ## The Armijo Goldstein line search scheme + accelerated gradient descent
  if (opts$mFlag == 0 && opts$lFlag == 0) {
    
    # this flag tests whether the gradient step only changes a little
    bFlag <- 0 
    
    # the intial guess of the Lipschitz continuous gradient
    L <- 1 / n + rsL2
    
    # the product between weight and y
    weighty <- weight * y
    
    #initial assignments
    betap <- beta
    xbetap <- xbeta
    betabetap <- numeric(p)
    cp <- c; ccp <- 0

    alphap <- 0; alpha <- 1
    
    for (iterStep in 1:opts$maxIter) {
      
      ## --------------------------- step 1 ---------------------------
      ##      compute search point s based on xp and x (with beta)
      bet <- (alphap - 1) / alpha
      s <- beta + bet * betabetap
      sc <- c + bet * ccp
      
      ## --------------------------- step 2 ---------------------------
      ##  line search for L and compute the new approximate solution x
      
      # compute xs = x * s
      xs <- xbeta + bet * (xbeta - xbetap)
      
      aa <- -y * (xs + sc)
      
      # fun.s is the logistic loss at the search point
      bb <- pmax(aa, 0)
      fun.s <- as.double( crossprod(weight, (log(exp(-bb) + exp(aa - bb)) + bb)) ) + 
        ( rsL2 / 2 ) * as.double(crossprod(s))
      
      # compute prob=[p_1;p_2;...;p_n]
      prob <- 1 / (1 + exp(aa))

      b <- -weighty * (1 - prob)

      #the gradient of c
      gc <- sum(b) 
      
      #compute g= xT b, the gradient of beta
      if (opts$nFlag == 0) {
        g <- crossprod(x, b)
      } else if (opts$nFlag == 1) {
        g <- (crossprod(x, b) - sum(b) * mu) / nu
      } else {
        invNu <- b / nu
        g <- crossprod(x, invNu) - sum(invNu) * mu
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
        
        
        if (is.null(groups)) {
          res <- flsa(v, z0, lambda / L, lambda2 / L, p,
                      1000, 1e-9, 1, 6)
          beta <- res[[1]]
          z0 <- res[[2]]
          infor <- res[[3]]
        } else {
          
          if (any(is.na(groups))) {
            ## don't apply fused lasso penalty
            ## to variables with group == NA 
            gr.idx <- which(is.na(groups))
            gr.p <- length(gr.idx)
            if (any(gr.idx == 1)) {
              gr.idx.z <- gr.idx[gr.idx != 1] - 1
            } else {
              gr.idx.z <- gr.idx[-gr.p]
            }
            
            res <- flsa(v[gr.idx], z0[gr.idx.z], lambda / L, 0, gr.p,
                        1000, 1e-9, 1, 6)
            
            beta[gr.idx] <- res[[1]]
            z0[gr.idx.z] <- res[[2]]
            infor <- res[[3]]
          }
          
          for (t in 1:length(unique.groups)) {
            gr.idx <- which(groups == unique.groups[t])
            gr.p <- length(gr.idx)
            if (any(gr.idx == 1)) {
              gr.idx.z <- gr.idx[gr.idx != 1] - 1
            } else {
              gr.idx.z <- gr.idx[-gr.p]
            }
            
            res <- flsa(v[gr.idx], z0[gr.idx.z], lambda / L, lambda2 / L, gr.p,
                        1000, 1e-9, 1, 6)
            
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
            
            beta[gr.idx] <- newbeta
            z0[gr.idx.z] <- res[[2]]
            infor <- res[[3]]
          }
        }
        
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
        
        aa <- -y * (xbeta + c)
        
        # fun.beta is the logistic loss at the new approximate solution
        bb <- pmax(aa, 0)
        
        fun.beta <- as.double( crossprod(weight, (log(exp(-bb) + exp(aa - bb)) + bb)) ) + 
          ( rsL2 / 2 ) * as.double(crossprod(beta))

        r.sum <- (as.double(crossprod(v)) / 2 + (c - sc)^2) / 2
        #l.sum <- fun.beta - fun.s - as.double(crossprod(v, g)) - (c - sc) * gc
        fzp.gamma <- fun.s + as.double(crossprod(v, g)) + (c - sc) * gc + L * r.sum
        
        
        if (r.sum <= 1e-18) {
          #this shows that the gradient step makes little improvement
          bFlag <- 1
          break
        }
        
        ## the condition is fun.beta <= fun.s + v'* g + c * gc
        ##                           + L/2 * (v'*v + (c-sc)^2 )
        
        #if (l.sum <= r.sum * L) {
        if (fun.beta <= fzp.gamma) {
          break
        } else {
          #L <- max(2 * L, l.sum / r.sum)
          L <- 2 * L
        }
        
      } # end while loop
      
      ## --------------------------- step 3 ---------------------------
      ##      update alpha and alphap, and check for convergence
      alphap <- alpha
      alpha <- (1 + sqrt(4 * alpha * alpha + 1)) / 2
      
      ## store values for L
      ValueL[iterStep] <- L
      
      betabetap <- beta - betap
      ccp <- c - cp
      
      # evaluate fused and group lasso 
      # penalty-terms
      if (!is.null(groups)) {
        fused.pen <- group.pen <- 0
        for (t in 1:length(unique.groups)) {
          gr.idx <- which(groups == unique.groups[t])
          gr.p <- length(gr.idx)
          if (gr.p > 1) {
            fused.pen <- fused.pen + sum(abs(beta[gr.idx[2:(gr.p)]] - beta[gr.idx[1:(gr.p - 1)]]))
            group.pen <- group.pen + sqrt(sum(beta[gr.idx] ^ 2) * gr.p)
          }
        }
        pens <- lambda2 * fused.pen + lambda.group * group.pen
      } else {
        pens <- lambda2 * sum(abs(beta[2:p] - beta[1:(p-1)]))
      }
      
      funVal[iterStep] <- fun.beta + lambda * sum(abs(beta)) + pens
      
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
  
  list(beta = drop(beta), intercept = c, funVal = funVal, ValueL = ValueL)
}