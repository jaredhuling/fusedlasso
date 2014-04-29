fusedlasso
====

__WARNING:__ this is far from a finished product (family == "gaussian" does not work properly yet)

an R translation of the code for the Efficient Fused Lasso Algorithm in 

_An Efficient Algorithm for a Class of Fused Lasso Problems_, Liu et al, 2010

http://www.public.asu.edu/~jye02/Publications/Papers/rp589f-liu.pdf

http://www.public.asu.edu/~jye02/Software/SLEP/index.htm


C code by the aforementioned authors (not me)


## Installation

**fusedlasso** is not available on CRAN, but can be installed using the R package **devtools**. **fusedlasso** can be installed with the following R code:

```r
devtools::install_github("jaredhuling/fusedlasso")
library(fusedlasso)
```

## Example

```r
nobs <- 10000
nvars <- 50

#generate data
set.seed(123)
true.beta <- rnorm(nvars) * rbinom(nvars, 1, 0.25)
true.beta[1:3] <- c(1, 1.06, 0.95)
x <- matrix(rnorm(nobs * nvars), ncol = nvars)

#generate binary outcome
log.p.ratio <- x %*% true.beta
prob.y.1 <- 1 / (1 + exp(-log.p.ratio))
y <- rbinom(nobs, 1, prob = prob.y.1)
y1 <- ifelse(y == 0, -1, y)

#fit fused lasso logistic model
res <- fusedlasso(x, y1, lambda.lasso = 0.005, lambda.fused = 0.01, family = "binomial")
round(res$beta, 5)

```