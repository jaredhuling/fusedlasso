

fusedLeastR <- function(x, y, lambda, lambda.group = 0, groups = NULL, opts=NULL) {
  
  # Function fusedLeastR
  #      Least Squares Loss with the Fused Lasso Penalty
  #
  ## Problem
  #
  #  min  1/2 || X b - y||^2  + lambda * ||b||_1 +   lambda_2 ||Rb||_1
  #                               ## + 1/2 rsL2 * ||b||_2^2
  #
  #  ##By default, rsL2=0.
  #  ##When rsL2 is nonzero, this correspons the well-know elastic net.
  #
  
  sz <- dim(x)
  n <- sz[1]
  p <- sz[2]

  stopifnot(lambda > 0 && lambda.group >= 0)
  
  # run sllOpts to set default values (flags)
  opts <- sllOpts(opts)
  if (lambda.group > 0) {
    opts$tol <- 1e-10
  }
  
  # if groups are given, get unique groups
  if (!is.null(groups)) {
    unique.groups <- sort(unique(groups[!is.na(groups)]))
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

  #### Starting point initialization
  
  ## compute X'y
  if (opts$nFlag == 0) {
    xTy <- crossprod(x, y)
  } else if (opts$nFlag == 1) {
    xTy <- (crossprod(x, y) - sum(y) * mu) / nu
  } else {
    invNu <- y / nu
    xTy <- crossprod(x, invNu) - sum(invNu) * mu
  }
  
  ## L1 norm regularization
  if (opts$rFlag == 0) {
    if (is.null(opts$fusedPenalty)) {
      lambda2 <- 0
    } else {
      lambda2 <- opts$fusedPenalty * n
    }
    lambda <- lambda * n
    lambda.group <- lambda.group * n
  } else {
    # lambda here is the scaling factor lying in [0,1]
    if (lambda < 0 || lambda > 1) {
      stop('opts$rFlag=1, and lambda should be in [0,1]')
    }
    
    lambda.max <- max(abs(xTy))
    lambda <- lambda * lambda.max
    lambda.group <- lambda.group * lambda.max
    
    if (is.null(opts$fusedPenalty)) {
      lambda2 <- 0
    } else {
      lambda2 <- lambda.max * opts$fusedPenalty
    }
    
  }
  
  ## initialize a starting point
  if (opts$init == 2) {
    b <- numeric(p)
  } else {
    if (!is.null(opts$x0)) {
      b <- opts$x0
      if (length(b) != p) {
        stop("Input x0 is of wrong length")
      }
    } else {
      b <- xTy
    }

  }
  
  
  ## compute x b
  if (opts$nFlag == 0) {
    xb <- x %*% b
  } else if (opts$nFlag == 1) {
    invNu <- b / nu
    mu.invNu <- as.double(crossprod(mu, invNu))
    xb <- x %*% invNu - rep(mu.invNu, n)
  } else {
    mub <- as.double(crossprod(mu, b))
    xb <- (x %*% b - rep(mub, n)) / nu
  }
  if (opts$init == 0) {
    ## If .init=0, we set x=ratio*x by "initFactor"
    
    b.norm <- sum(abs(b))
    b.2norm <- as.double(crossprod(b))
    
    if (b.norm >= 1e-6) {
      ratio <- initFactor(b.norm, xb, y, lambda, "LeastR", 0, b.2norm)
      b <- ratio * b
      xb <- ratio * xb
    }
  }
  
  if (is.null(opts$rStartNum)) {
    opts$rStartNum <- opts$maxIter
  } else {
    if (opts$rStartNum <=0) {
      opts$rStartNum <- opts$maxIter
    }
  }
  
  
  #### MAIN CODE ####
  
  ## z0 is the starting point for flsa
  z0 <- numeric((p-1))
  
  ## this flag tests whether the gradient 
  ## step only changes a little
  bFlag <- 0
  
  ## L=1 + rsL2;
  L <- 1
  
  ValueL <- funVal <- numeric(opts$maxIter)
  
  # We assume that the maximum eigenvalue of A'A is over 1
  
  ## The Armijo Goldstein line search scheme + accelerated gradient descent
  if (opts$lFlag == 0) {
    
    ## assign bp with b, and xbp with xb
    bp <- b; xbp <- xb; bbp <- numeric(p)
        
    ## alphap and alpha are used for computing the weight in forming search point
    alphap <- 0; alpha <- 1;
    
    for (iterStep in 1:opts$maxIter) {
      ## --------------------------- step 1 ---------------------------
      ##      compute search point s based on xp and x (with beta)
      beta <- (alphap - 1) / alpha
      s <- b + beta * bbp
      
      ## --------------------------- step 2 ---------------------------
      ##  line search for L and compute the new approximate solution x
      
      ## compute the gradient (g) at s
      xs <- xb + beta * (xb - xbp)
      
      ## compute xT, xs
      if (opts$nFlag == 0) {
        xTxs <- crossprod(x, xs)
      } else if (opts$nFlag == 1) {
        xTxs <- (crossprod(x, xs) - sum(xs) * mu) / nu
      } else {
        invNu <- xs / nu
        xTxs <- crossprod(x, invNu) - sum(invNu) * mu
      }
      
      ## obtain the gradient g
      g <- (xTxs - xTy)
      
      ## copy b and xb to bp and xbp
      bp <- b; xbp <- xb
      
      while (TRUE) {
        ## let's walk in a step in the antigradient of s to get v
        ## and then do the l1-norm regularized projection
        v <- s - g / L
        
        if (is.null(groups)) {
          res <- flsa(v, z0, lambda / L, lambda2 / L, p,
                      1000, 1e-8, 1, 6)
          b <- res[[1]]
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
                        1000, 1e-8, 1, 6)
            b[gr.idx] <- res[[1]]
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
                newbeta = (pmax(nm - lambda.group * sqrt(gr.p), 0) / nm) * res[[1]]
              }
              end
            } else {
              newbeta <- res[[1]]
            }
            
            b[gr.idx] <- newbeta
            z0[gr.idx.z] <- res[[2]]
            infor <- res[[3]]
          }
        }
        
        ## the difference between the new approximate 
        ## solution x and the search point s
        v <- b - s
        
        ## compute x b
        if (opts$nFlag == 0) {
          xb <- x %*% b
        } else if (opts$nFlag == 1) {
          invNu <- b / nu
          mu.invNu <- as.double(crossprod(mu, invNu))
          xb <- x %*% invNu - rep(mu.invNu, n)
        } else {
          mub <- as.double(crossprod(mu, b))
          xb <- (x %*% b - rep(mub, n)) / nu
        }
        
        xv <- xb - xs
        
        r.sum <- as.double(crossprod(v))
        l.sum <- as.double(crossprod(xv))
        
        if (r.sum < 1e-18) {
          # this shows that the gradient step makes little improvement
          bFlag <- 1
          break
        }
        
        ## the condition is ||Av||_2^2 <= (L - rsL2) * ||v||_2^2
        if (l.sum <= r.sum * L) {
          break
        } else {
          L <- max(2 * L, l.sum / r.sum)
        }
        
      } # end while loop
      
      ValueL[iterStep] <- infor[2]
      
      ## --------------------------- step 3 ---------------------------
      ##      update alpha and alphap, and check for convergence
      alphap <- alpha
      alpha <- (1 + sqrt(4 * alpha * alpha + 1)) / 2
      
      bbp <- b - bp
      xby <- xb - y
      
      # evaluate fused and group lasso 
      # penalty-terms
      if (!is.null(groups)) {
        fused.pen <- group.pen <- 0
        for (t in 1:length(unique.groups)) {
          gr.idx <- which(groups == unique.groups[t])
          gr.p <- length(gr.idx)
          if (gr.p > 1) {
            fused.pen <- fused.pen + sum(abs(b[gr.idx[2:(gr.p)]] - b[gr.idx[1:(gr.p - 1)]]))
            group.pen <- group.pen + sqrt(sum(b[gr.idx] ^ 2) * gr.p)
          }
        }
        pens <- lambda2 * fused.pen + lambda.group * group.pen
      } else {
        pens <- lambda2 * sum(abs(b[2:p] - b[1:(p-1)]))
      }
      
      funVal[iterStep] <- as.double(crossprod(xby)) / (2) + lambda * sum(abs(b)) + pens


      
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
        bp <- b; xbp <- xb; L <- 1
        bbp <- numeric(p)
      }
      
    } #end loop
    funVal <- funVal[1:iterStep]
    ValueL <- ValueL[1:iterStep]
    
  } # end accelerated gradient if-statement
  
  
  ## the line search scheme by Nesterov
  if (opts$lFlag == 1) {
    
    ## compute xTxb
    if (opts$nFlag == 0) {
      xTxb <- crossprod(x, xb)
    } else if (opts$nFlag == 1) {
      xTxb <- (crossprod(x, xb) - sum(xb) * mu) / nu
    } else {
      invNu <- ab / nu
      xTxb <- crossprod(x, invNu) - sum(invNu) * mu
    }
    g.b <- xTxb - xTy
    
    ## initialization
    v <- b      #v is the minimizer of \psi_i(x)
    g.v <- g.b  #g.v is the gradient at v
    sum.g <- b  #sum_g contains the summation of b0 - \sum_{i >=1} g_{b_i}
    sum.a <- 0  #a is a nonnegative tuning parameter, 
                #and sum.a is the summation of all the
                #a's during the iteration  
    
    for (iterStep in 1:opts$maxIter) {
      while (TRUE) {
        # compute a, which is the solution to the quadratic function
        a <- (1 + sqrt(1 + 2 * L * sum.a)) / L
        
        # compute the search point s
        s <- (sum.a / (sum.a + a)) * b + (a / (sum.a + a)) * v
        
        # the gradient at the search point s
        g.s <- (sum.a / (sum.a + a)) * g.b + (a / (sum.a + a)) * g.v
        
        # compute the new solution bnew
        u <- s - g.s / L      # u is a gradient based on s
        
        #obtain bnew by projection
        z0 <- z0 / L
        
        res <- flsa(u, z0, lambda / L, lambda2 / L, p,
                    1000, 1e-8, 1, 6)
        
        bnew <- res[[1]]
        z <- res[[2]]
        infor <- res[[3]]
        z0 <- z * L
        
        # compute x bnew
        if (opts$nFlag == 0) {
          xbnew <- x %*% bnew
        } else if (opts$nFlag == 1) {
          invNu <- bnew / nu
          mu.invNu <- as.double(crossprod(mu * invNu))
          xbnew <- x %*% invNu - rep(mu.invNu, n)
        } else {
          mub <- as.double(crossprod(mu, b))
          xbnew <- (x %*% bnew - rep(mub, n)) / nu
        } 
        
        xb <- xbnew
        # compute xT xb
        if (opts$nFlag == 0) {
          xTxb <- crossprod(x, xb)
        } else if (opts$nFlag == 1) {
          xTxb <- (crossprod(x, xb) - sum(xb) * mu) / nu
        } else {
          invNu <- xb / nu
          xTxb <- crossprod(x, invNu) - sum(invNu) * mu
        } 
        
        #g.bnew is the gradient at bnew
        g.bnew <- xTxb - xTy
        
        # test whether L is appropriate
        # cG denotes the composite gradient
        cG <- L * (s - bnew) + g.bnew - g.s
        val.left <- as.double(crossprod(cG, (s - bnew)))
        val.right <- as.double(crossprod(cG))

        if (crossprod(g.bnew) < 1e-18) {
          # this shows that the gradient step makes little improvement
          bFlag <- 1
          break
        }
        
        if (val.left * L >= val.right) {
          ValueL[iterStep] <- L
          L <- L / 1.1
          break
        } else {
          L <- max(2 * L)
        }
        
      }
      
      # compute v as the minizer of \psi_k (b)
      sum.g <- sum.g - g.bnew * a
      sum.a <- sum.a + a
      
      # obtain v by projection
      z0 <- z0 * sum.a
      
      res <- flsa(sum.g, z0, lambda * sum.a, lambda2 * sum.a, p,
                  1000, 1e-8, 1, 6)
      
      v <- res[[1]]
      z <- res[[2]]
      infor <- res[[3]]
      z0 <- z / sum.a
      
      # compute xv
      if (opts$nFlag == 0) {
        xv <- x %*% v
      } else if (opts$nFlag == 1) {
        #######
        ####### This doesn't work yet
        #######
        invNu <- v / nu
        mu.invNu <- mu * invNu
        xv <- x %*% invNu - rep(mu.invNu, each = n)
      } else {
        #######
        ####### This doesn't work yet
        #######
        xv <- (x %*% v - rep(mu * v, each = m)) / nu
      }
      
      # compute xT xv
      if (opts$nFlag == 0) {
        g.v <- crossprod(x, xv)
      } else if (opts$nFlag == 1) {
        g.v <- (crossprod(x, xv) - sum(xv) * mu) / nu
      } else {
        invNu <- xv / nu
        g.v <- crossprod(x, invNu) - sum(invNu) * mu
      }
      g.v <- g.v - xTy 
      
      bbp <- bnew - b
      norm.bp <- norm(b, type = "F")
      norm.bbp <- norm(bbp, type = "F")
      b <- bnew
      g.b <- g.bnew
      
      # Compute the objective function value
      xby <- xb - y
      
      funVal[iterStep] <- as.double(crossprod(xby)) / 2 + 
        sum(abs(b)) * lambda + lambda2 * sum(abs( b[2:p] - b[1:(p-1)] ))
      
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
        v <- b
        g.v <- g.b
        sum.g <- b
        sum.a <- 0
      }
      
      
    } #end iterStep loop
    funVal <- funVal[1:iterStep]
    ValueL <- ValueL[1:iterStep]
  } # end Nesterov line-search
  
  
  
  list(beta = drop(b), funVal = funVal, ValueL = ValueL)
}





