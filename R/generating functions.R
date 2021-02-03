### Normal ---------------------------------------------------------------------
NO <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu = NULL){
    stats::dnorm(sqrt(u))
  }

  ## Cumulative distribution function
  out$p <- function(q, nu = NULL){
    stats::pnorm(q)
  }

  ## Quantile function
  out$q <- function(p, nu = NULL){
    stats::qnorm(p)
  }

  ## Extra parameter indicator
  out$extraP <- FALSE

  ## Wheighting function
  out$weigh <- function(z, nu = NULL) 1

  out
}

### Student-t ------------------------------------------------------------------
ST <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu){
    stats::dt(sqrt(u), nu)
  }

  ## Cumulative distribution function
  out$p <- function(q, nu){
    stats::pt(q, nu)
  }

  ## Quantile function
  out$q <- function(p, nu){
    stats::qt(p, nu)
  }

  ## Extra parameter indicator
  out$extraP <- TRUE

  ## Wheighting function
  out$weigh <- function(z, nu) (nu + 1)/(nu + z^2)

  out
}

### Power exponential ----------------------------------------------------------
# From gamlss.dist package (dPE, pPE, qPE)
PE <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu){
    x <- sqrt(u)

    if (any(nu < 0))
      stop("nu must be positive\n")

    log.c <- 0.5 * (-(2/nu) * log(2) + lgamma(1/nu) - lgamma(3/nu))
    c <- exp(log.c)

    log.lik <-  log(nu) - log.c - (0.5 * (abs(x/c)^nu)) -
      (1 + (1/nu)) * log(2) - lgamma(1/nu)

    return(exp(log.lik))

  }

  ## Cumulative distribution function
  out$p <- function(q, nu){
    if (any(nu < 0))
      stop("nu must be positive\n")

    log.c <- 0.5 * (-(2/nu) * log(2) + lgamma(1/nu) - lgamma(3/nu))
    c <- exp(log.c)
    s <- 0.5 * ((abs(q/c))^nu)

    cdf <- 0.5 * (1 + stats::pgamma(s, shape = 1/nu, scale = 1) * sign(q))

    id <- which(nu > 10000, arr.ind = TRUE)
    cdf[id] <- (q[id] + sqrt(3))/sqrt(12)

    cdf
  }

  ## Quantile function
  out$q <- function(p, nu){
    lower.tail <- TRUE
    log.p <- FALSE

    if (any(nu < 0))
      stop("nu must be positive\n")

    if (any(p < 0) | any(p > 1))
      stop("p must be between 0 and 1\n")

    log.c <- 0.5 * (-(2/nu) * log(2) + lgamma(1/nu) - lgamma(3/nu))
    c <- exp(log.c)
    suppressWarnings(s <- stats::qgamma((2 * p - 1) * sign(p - 0.5),
                                        shape = (1/nu), scale = 1))
    z <- sign(p - 0.5) * ((2 * s)^(1/nu)) * c
    z
  }

  ## Extra parameter indicator
  out$extraP <- TRUE

  ## Wheighting function
  out$weigh <- function(z, nu){
    log.c <- 0.5 * (-(2/nu) * log(2) + lgamma(1/nu) - lgamma(3/nu))
    c <- exp(log.c)
    (nu * (z^2)^(nu/2 - 1)) / (2 * c^nu)
  }

  out
}

### Cauchy ---------------------------------------------------------------------
CA <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu = NULL){
    stats::dcauchy(sqrt(u))
  }

  ## Cumulative distribution function
  out$p <- function(q, nu = NULL){
    stats::pcauchy(q)
  }

  ## Quantile function
  out$q <- function(p, nu = NULL){
    stats::qcauchy(p)
  }

  ## Extra parameter indicator
  out$extraP <- FALSE

  ## Wheighting function
  out$weigh <- function(z, nu = NULL) 2 / (1 + z^2)

  out
}

### Double exponential ---------------------------------------------------------
# From rmulti package
DE <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu = NULL){
    (sqrt(2) / 2) * exp(-sqrt(2 * u))
  }

  ## Cumulative distribution function
  out$p <- function(q, nu = NULL){
    u <- sqrt(2) * q
    t <- exp(-abs(u)) / 2
    ifelse(u < 0, t, 1 - t)
  }

  ## Quantile function
  out$q <- function(p, nu = NULL){
    h <- function(y){
      u <- sqrt(2) * y
      t <- exp(-abs(u)) / 2
      ifelse(u<0, t, 1 - t) - p
    }
    if(any(p < 0 | p > 1))
      stop("p must lie between 0 and 1\n")

    ifelse(p < 0.5, 1 / sqrt(2) * log(2 * p),-1 / sqrt(2) * log(2 * (1 - p)))
  }

  ## Extra parameter indicator
  out$extraP <- FALSE

  ## Wheighting function
  out$weigh <- function(z, nu = NULL) sqrt(2) / abs(z)

  out
}

### Logistic -------------------------------------------------------------------
LO <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu = NULL){
    stats::dlogis(sqrt(u))
  }

  ## Cumulative distribution function
  out$p <- function(q, nu = NULL){
    stats::plogis(q)
  }

  ## Quantile function
  out$q <- function(p, nu = NULL){
    if(any(p < 0 | p > 1))
      stop("p must lie between 0 and 1\n")

    stats::qlogis(p)
  }

  ## Extra parameter indicator
  out$extraP <- FALSE

  ## Wheighting function
  out$weigh <- function(z, nu = NULL){
    (exp(-abs(z)) - 1) / (abs(z) * (exp(-abs(z)) + 1))
  }

  out
}
