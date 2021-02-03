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

### Canonical slash ------------------------------------------------------------
CSL <- function(){
  out <- list()

  ## Generating function
  out$d <- function(u, nu = NULL){
    p0 <- rep(1 / (2 * sqrt(2 * pi)), length(u[u == 0]) )
    p <- (1 - exp(-u[u > 0]/2)) / (sqrt(2 * pi) * u[u > 0])

    pdf <- c(p0, p)
    index <- c(which(u == 0), which(u > 0))

    pdf[sort(index, index.return = T)$ix]
  }

  ## Cumulative distribution function
  out$p <- function(q, nu = NULL){
    p0 <- rep(0.5, length(q[q == 0]) )
    p <- pnorm(q[q != 0]) - (dnorm(0) - dnorm(q[q != 0])) / q[q != 0]

    cdf <- c(p0, p)
    index <- c(which(q == 0), which(q != 0))

    cdf[sort(index, index.return = T)$ix]
  }

  ## Quantile function
  out$q <- function(p, nu = NULL){
    if(any(p < 0 | p > 1))
      stop("p must lie between 0 and 1\n")

    qtf <- function(p){
      if (length(p) >= 1){
        obj <- function(q){
          pnorm(q) - (dnorm(0) - dnorm(q)) /
            q - p
        }

        #jac <- function(q){
        #  dnorm(0)/(q^2) - dnorm(q)
        #}

        nleqslv::nleqslv(2 * qnorm(p), obj)$x
      }else{
        numeric(0)
      }
    }

    q0 <- rep(0, length(p[p == 0.5]) )
    q <- as.numeric(apply(matrix(p[p != 0.5], ncol = 1), 1, qtf))

    qtf <- c(q0, q)
    index <- c(which(p == 0.5), which(p != 0.5))

    qtf[sort(index, index.return = T)$ix]
  }

  ## Extra parameter indicator
  out$extraP <- FALSE

  ## Wheighting function
  out$weigh <- function(z, nu = NULL){
    2/(z^2) - exp(-z^2 / 2) / (1 - exp(- z^2 / 2))
  }

  out
}

### Slash ----------------------------------------------------------------------
SL <- function(){
  out <- list()

  ## Incomplete gamma function
  ig <- function(a, x) gamma(a) * pgamma(x, a)

  ## Generating function
  out$d <- function(u, nu){

    if (any(nu < 0))
      stop("nu must be positive\n")

    if (length(nu) == 1) nu <- rep(nu, length(u))

    p0 <- rep(nu[u == 0] / ((nu[u == 0] + 1) * sqrt(2 * pi)), length(u[u == 0]) )
    p <- ig((nu[u > 0] + 1)/2, u[u > 0]/2) * nu[u > 0] * 2^(nu[u > 0]/2 - 1) /
      (sqrt(pi) * (u[u > 0]^((nu[u > 0] + 1)/2)))

    pdf <- c(p0, p)
    index <- c(which(u == 0), which(u > 0))

    pdf[sort(index, index.return = T)$ix]
  }

  ## Cumulative distribution function
  out$p <- function(q, nu){

    if (any(nu < 0))
      stop("nu must be positive\n")

    if (length(nu) == 1) nu <- rep(nu, length(q))

    p0 <- rep(0.5, length(q[q == 0]))

    CDF <- function(input){
      q <- input[1]
      nu <- input[2]

      if (length(q) >= 1){

        integrand <- function(x){
          nu * x^(nu-1) * pnorm(x * q)
        }

        integrate(integrand, 0, 1)$val
      }else{
        numeric(0)
      }
    }

    p <- as.numeric(apply(matrix(c(q[q != 0], nu[q != 0]), ncol = 2), 1, CDF))

    cdf <- c(p0, p)
    index <- c(which(q == 0), which(q != 0))

    cdf[sort(index, index.return = T)$ix]
  }

  ## Quantile function
  out$q <- function(p, nu){
    if(any(p < 0 | p > 1))
      stop("p must lie between 0 and 1\n")

    if (any(nu < 0))
      stop("nu must be positive\n")

    if (length(nu) == 1) nu <- rep(nu, length(p))

    qtf <- function(input){
      p <- input[1]
      nu <- input[2]

      if (!is.na(p)){
        obj <- function(q){
          out$p(q, nu) - p
        }

        nleqslv::nleqslv(qnorm(p) / (0.5^(1/nu)), obj)$x
      }else{
        numeric(0)
      }
    }

    q0 <- rep(0, length(p[p == 0.5]) )
    q <- as.numeric(apply(matrix(c(p[p != 0.5], nu[p != 0.5]), ncol = 2), 1, qtf))

    qtf <- c(q0, q)
    index <- c(which(p == 0.5), which(p != 0.5))

    qtf[sort(index, index.return = T)$ix]
  }

  ## Extra parameter indicator
  out$extraP <- TRUE

  ## Wheighting function
  out$weigh <- function(z, nu){
    2 * ig( (nu + 3)/2, z^2 / 2 ) / ( z^2 * ig( (nu + 1)/2, z^2 / 2 ))
  }

  out
}

