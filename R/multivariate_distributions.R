### Multivariate normal --------------------------------------------------------
# From 'mvtnorm' package
gaussian <- function(){
  out <- list()

  ### Multivariate density
  out$Md <- function(x, P = diag(d), checkSymmetry = TRUE){

    if (is.vector(x))
      x <- matrix(x, ncol = length(x))

    d <- ncol(x)

    if (!missing(P)) {
      if (d != ncol(P))
        stop("\nx and P have non-conforming size\n")

      if (checkSymmetry && !isSymmetric(P, tol = sqrt(.Machine$double.eps),
                                        check.attributes = FALSE))
        stop("\nP must be a symmetric matrix\n")

    }

    dec <- tryCatch(chol(P), error = function(e) e)
    if (inherits(dec, "error")) {
      x.is.mu <- colSums(t(x) != rep(0, d)) == 0
      logretval <- rep.int(-Inf, nrow(x))
      logretval[x.is.mu] <- Inf
    }else {
      tmp <- backsolve(dec, t(x), transpose = TRUE)
      rss <- colSums(tmp^2)
      logretval <- -sum(log(diag(dec))) - 0.5 * d * log(2 * pi) - 0.5 * rss
    }

    names(logretval) <- rownames(x)
    exp(logretval)
  }

  ### Random generation
  out$Mr <- function(n, P = diag(d), checkSymmetry = TRUE){

    d <- dim(P)[1]

    if (checkSymmetry && !isSymmetric(P, tol = sqrt(.Machine$double.eps),
                                      check.attributes = FALSE)) {
      stop("\nP must be a symmetric matrix\n")
    }

    ev <- eigen(P, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("\nP is numerically not positive semidefinite\n")
    }
    R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))


    retval <- matrix(stats::rnorm(n * ncol(P)), nrow = n) %*% R
    retval <- sweep(retval, 2, rep(0, d), "+")
    retval
  }

  ## Name
  out$name <- "Gaussian"

  out

}

### Multivariate Student-t -----------------------------------------------------
# From 'mvtnorm' package
st <- function(){
  out <- list()

  ### Multivariate density
  out$Md <- function(x, P = diag(d), df = 4, checkSymmetry = TRUE){

    if (is.vector(x))
      x <- matrix(x, ncol = length(x))

    d <- ncol(x)

    if (!missing(P)) {
      if (d != ncol(P))
        stop("\nx and P have non-conforming size\n")

      if (checkSymmetry && !isSymmetric(P, tol = sqrt(.Machine$double.eps),
                                        check.attributes = FALSE))
        stop("\nP must be a symmetric matrix\n")

    }

    dec <- tryCatch(chol(P), error = function(e) e)
    if (inherits(dec, "error")) {
      x.is.mu <- colSums(t(x) != rep(0, d)) == 0
      logretval <- rep.int(-Inf, nrow(x))
      logretval[x.is.mu] <- Inf
    }else {
      R.x_m <- backsolve(dec, t(x), transpose = TRUE)
      rss <- colSums(R.x_m^2)
      logretval <- lgamma((d + df)/2) - (lgamma(df/2) + sum(log(diag(dec))) +
                                           d/2 * log(pi * df)) - 0.5 * (df + d) * log1p(rss/df)
    }

    names(logretval) <- rownames(x)
    exp(logretval)
  }

  ### Random generation
  out$Mr <- function(n, P = diag(d), df = 4, checkSymmetry = TRUE){

    d <- dim(P)[1]

    if (checkSymmetry && !isSymmetric(P, tol = sqrt(.Machine$double.eps),
                                      check.attributes = FALSE)) {
      stop("\nP must be a symmetric matrix\n")
    }

    ev <- eigen(P, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("\nP is numerically not positive semidefinite\n")
    }

    sims <- gaussian()$Mr(n, P = P)/sqrt(stats::rchisq(n, df)/df)
    sweep(sims, 2, rep(0, d), "+")
  }

  ## Name
  out$name <- "Student's t"

  out

}
