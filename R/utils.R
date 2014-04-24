
value.func <- function(A, predicted, y) {
  A.char <- levels(A)[A]
  ind <- 1 * (A.char == predicted)
  sum(3 * y * ind) / length(A)
}


logLikNull <- function(y, A, weights) {
  require(mlogit)
  mdat <- data.frame(A = A, weights = weights)
  rownames(mdat) <- 1:length(weights)
  mdat2 <- mlogit.data(mdat, choice = "A", shape = "wide")
  mfit <- mlogit(A ~ 1, weights = weights, data = mdat2)
  logLik(mfit)[1]
}

logLikGlmnet <- function(y, obj, A, weights) {
  logLikSat <- (obj$nulldev / 2) + logLikNull(y, A, weights)
  logLikSat - deviance(obj) / 2
}

aic <- function(obj, y, A, weights) {
  stopifnot(inherits(obj, "glmnet"))
  -2 * logLikGlmnet(y, obj, A, weights) + 2 * obj$df
}

df.msgl <- function(obj) {
  betas <- array(0, dim = length(obj$beta))
  for (i in 1:length(obj$beta)) {
    betas[i] <- sum(colSums(as(obj$beta[[i]], "matrix")) != 0)
  }
  betas
}

aic.msgl <- function(obj) {
  2 * (df.msgl(obj) + obj$loss * length(obj$classes.true))
}

bic <- function(obj, y, A, weights) {
  stopifnot(inherits(obj, "glmnet"))
  -2 * logLikGlmnet(y, obj, A, weights) + obj$df * (log(obj$nobs) + log(2 * pi))
}


predTreatment <- function(obj, x) {
  if (inherits(obj, "cv.glmnet")) {
    obj <- obj$glmnet.fit
  }
  coef.vals <- array(0, dim = c(length(obj$beta), nrow(x), ncol(obj$beta[[1]])))
  for (i in 1:length(obj$beta)) {
    coef.vals[i, , ] <- drop(cbind(rep(1, nrow(x)), x) %*% rbind(obj$a0[i,], as(obj$beta[[i]], "matrix")))
  }
  apply(coef.vals, 2:3, function(x) as.character(which.max(x)))
}

modelMatrix <- function(data, interaction = FALSE, zero.one = TRUE) {
  
  varnames <- vn <- colnames(data)
  for (i in 1:ncol(data)){
    if (is.factor(data[[i]])) {
      if (length(levels(data[[i]])) < 3) {
        varnames[i] <- NA
      }
    } else {
      varnames[i] <- NA
    }
  }
  factornames <- varnames[!is.na(varnames)]
  
  x <- model.matrix(~ ., data = data)[,-1]
  p <- ncol(x)
  if (!zero.one) {
    x[x == 0] <- -1
  }
  x.ret <- array(0, dim = c(nrow(x), if(interaction)
  {choose(p, 2) + p} else {p}))
  colnames(x.ret) <- 1:ncol(x.ret)
  colx <- colnames(x)
  colnames(x.ret)[1:p] <- colx
  x.ret[,1:p] <- x
  groups <- rep(NA, ncol(x.ret))
  k <- 0
  factor.int.groups <- NULL
  for (i in 1:length(factornames)) {
    pattern <- paste(factornames[i], "[^:]*$", sep = "")
    factor.locs <- grep(pattern, colnames(x.ret))
    groups[factor.locs] <- i
  }
  if (interaction) {
    ig <- 1
    n.bad.ints <- 0
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        v1.root <- vn[which(pmatch(vn, colx[i]) == 1)]
        v2.root <- vn[which(pmatch(vn, colx[j]) == 1)]
        if (v1.root != v2.root) {
          k <- k + 1
          x.ret[,p + k] <- x[,i] * x[,j]
    
          if (any(!is.na(match(c(v1.root, v2.root), factornames)))) {
            gr <- paste(v1.root, v2.root, sep = "_")
            if (is.null(factor.int.groups)) {
              factor.int.groups[1] <- gr
            } else {
              if (is.na(match(gr, factor.int.groups))) {
                ig <- ig + 1
                factor.int.groups[ig] <- gr
              }
            }
            groups[p + k] <- match(gr, factor.int.groups) + length(factornames)
          }
          colnames(x.ret)[p + k] <- paste(colx[i], colx[j], sep = ":")
        } else {n.bad.ints <- n.bad.ints + 1}
      }
    }
    if (n.bad.ints > 0) {
      x.ret <- x.ret[,-((ncol(x.ret) - n.bad.ints +1):(ncol(x.ret)))]
      groups <- groups[-((length(groups) - n.bad.ints +1):(length(groups)))]
    }
  }

  names(groups) <- colnames(x.ret)
  attr(x.ret, "groups") <- groups
  x.ret
}

patientEffectData2d3 <- function(obj, x, patient.ind, lam.ind, 
                                 patient.names = NULL, json = TRUE) {
  require(rjson)
  n.dat <- length(patient.ind)
  ret.dat <- vector(length = n.dat, mode = "list")
  
  for (l in 1:n.dat) {
    effects <- trts <- NULL
    
    if (inherits(obj$model, "glmnet")) {
      ret <- vector(length = length(obj$model$beta), mode = "list")
      names(ret) <- names(obj$model$beta)
      for (i in 1:length(obj$model$beta)) {
        c.betas <- obj$model$beta[[i]][,lam.ind]
        nz.ind <- which(c.betas != 0)
        c.betas.nz <- c.betas[nz.ind]
        col.ind <- match(names(c.betas.nz), colnames(x))
        effects.tmp <- c(obj$model$a0[i, lam.ind], x[patient.ind[l], col.ind] * c.betas.nz)
        effects.tmp <- effects.tmp[order(abs(effects.tmp), decreasing = TRUE)]
        ret[[i]] <- effects.tmp
      }
    } else {
      full.beta <- as(obj$model$beta[[lam.ind]], "matrix")
      ret <- vector(length = nrow(obj$model$beta[[1]]), mode = "list")
      trt.names <- dimnames(full.beta)[[1]]
      for (i in 1:nrow(obj$model$beta[[1]])) {
        c.betas <- full.beta[i,-1]
        nz.ind <- which(c.betas != 0)
        c.betas.nz <- c.betas[nz.ind]
        col.ind <- match(names(c.betas.nz), colnames(x))
        effects.tmp <- c(full.beta[i,1], x[patient.ind[l], col.ind] * c.betas.nz)
        effects.tmp <- effects.tmp[order(abs(effects.tmp), decreasing = TRUE)]
        ret[[i]] <- effects.tmp
      }
    }
    
    ret.dat[[l]] <- ret
  }
  names(ret.dat) <- if (is.null(patient.names)) {
    paste("Patient", patient.ind, sep = "")
  } else {
    patient.names
  }
  if(json) {toJSON(ret.dat)} else {ret.dat}
}
