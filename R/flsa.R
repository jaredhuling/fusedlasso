
flsa <- function(x, z, 
                 lambda1, lambda2, n, 
                 maxStep, tol, tau, flag) {
  
  infor <- matrix(as.double(rep(0, 4)), ncol=1)
  res <- .C('R_flsa', 
            x=as.matrix(as.double(x), ncol = 1),
            z=as.matrix(as.double(z), ncol = 1), 
            infor=infor,
            as.double(lambda1), as.double(lambda2),
            as.integer(n), as.integer(maxStep),
            as.double(tol), as.integer(tau), 
            as.integer(flag),
            PACKAGE = "fusedlasso") 
  res
}