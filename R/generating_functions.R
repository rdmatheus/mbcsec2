#' @name BCSgen
#'
#' @title Box-Cox Symmetric Distribution Class Generators
#'
#' @description Current available generating distributions that can be used to
#'     fit a model belonging to the Box-Cox symmetric class of distributions.
#'
#' @param dist Character specification of the generating distribution, see
#'     details.
#' @param x A \code{"BCSgen"} object.
#' @param ... Further arguments for other specific methods.
#'
#' @details There are some distributions available for the density generating
#'     function. The following table display their names and their abbreviations
#'     to be passed to \code{BCSgen()}.
#'
#'   \tabular{lll}{
#'  \bold{Distribution}  \tab \bold{Abbreviation} \tab \bold{Does it have an extra parameter?}\cr
#'  Cauchy  \tab \code{"CA"}      \tab  no  \cr
#'  Canonical slash  \tab \code{"CSL"}      \tab  no  \cr
#'  Double exponential (Laplace)  \tab \code{"DE"}      \tab  no  \cr
#'  Logistic  \tab \code{"LO"}      \tab  no  \cr
#'  Normal  \tab \code{"NO"}      \tab  no  \cr
#'  Power exponential  \tab \code{"PE"}      \tab  yes  \cr
#'  Slash  \tab \code{"SL"}      \tab  yes  \cr
#'  Student-t  \tab \code{"ST"}      \tab  yes  \cr
#'  }
#'
#'
#' @return The function \code{BCSgen()} returns a list whose components are set
#'     of functions referring to the generating function of the Box-Cox
#'     symmetric class of distributions. More specifically, returns an
#'     \code{"BCSgen"} object with the following elements:
#'  \itemize{
#'    \item{d:}{ Density generating function \code{function(u, nu)}.}
#'    \item{p:}{ Cumulative distribution function \code{function(q, nu)}.}
#'    \item{q:}{ Quantile function \code{function(p, nu)}.}
#'    \item{extraP:}{ Logical; it receives \code{TRUE} if the generating
#'      distribution has an extra parameter.}
#'    \item{weigh:}{ Weighting function \code{function(z, nu)}.}
#'    \item{name:}{ Name of the distribution.}
#'  }
#'
#'  The arguments of the returned functions are:
#'  \describe{
#'    \item{\code{u, q}}{ Vector of positive quantiles.}
#'    \item{\code{p}}{ Vector of probabilities.}
#'    \item{\code{z}}{ A vector of transformed responses with the generalized
#'     Box-Cox transformation.}
#'    \item{\code{nu}}{ Extra parameter of the generating distribution, if it
#'     does not exist, it receives \code{NULL}.}
#'  }
#'
#' @references Ferrari, S. L., & Fumes, G. (2017). Box–Cox symmetric
#'  distributions and applications to nutritional data. AStA Advances in
#'  Statistical Analysis, 101(3), 321-344.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' \dontrun{
#' BCSgen("NO")
#'
#' curve(BCSgen("NO")$d(x^2), xlim = c(-5, 5),
#'       ylim = c(0, 0.7), ylab = "Density")
#' curve(BCSgen("CA")$d(x^2), add = TRUE, col = 2)
#' curve(BCSgen("CSL")$d(x^2), add = TRUE, col = 3)
#' curve(BCSgen("DE")$d(x^2), add = TRUE, col = 4)
#' curve(BCSgen("LO")$d(x^2), add = TRUE, col = 6)
#' legend("topright", c("Cauchy", "Canonical Slash",
#'                      "Laplace", "Logistic", "Normal"),
#'         col = c(2:4, 6, 1), lty = 1, bty = "n")
#' }
#'
#' @export
#'
BCSgen <- function(dist){
  dist <- match.fun(dist)
  dist <- eval(dist())
  class(dist) <- "BCSgen"
  dist
}

#' @rdname BCSgen
#' @export
print.BCSgen <- function(x, ...){

  gens <- c("CA", "CSL", "DE", "LO", "NO", "PE", "SL", "ST")
  GENs <- c("Cauchy", "Canonial slash", "Double exponential", "Logistic",
            "Normal", "Power exponential", "Slash", "Student-t")


  cat("----------------------------------------",
      "\nBox-Cox Symmetrical Generating Function",
      "\n----------------------------------------",
      "\nGenerating distribution:", x$name,
      "\nAbreviation:", gens[GENs == x$name],
      "\nExtra parameter:", ifelse(x$extraP, "yes", "no"),
      "\n---")
}

################################################################################
# Distributions:                                                               #
###############################################################################

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

  ## Name
  out$name <- "Normal"

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

  ## Name
  out$name <- "Student-t"

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

  ## Name
  out$name <- "Power exponential"

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

  ## Name
  out$name <- "Cauchy"

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

  ## Name
  out$name <- "Double exponential"

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

  ## Name
  out$name <- "Logistic"

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
    p <- stats::pnorm(q[q != 0]) - (stats::dnorm(0) - stats::dnorm(q[q != 0])) / q[q != 0]

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
          stats::pnorm(q) - (stats::dnorm(0) - stats::dnorm(q)) /
            q - p
        }

        #jac <- function(q){
        #  stats::dnorm(0)/(q^2) - stats::dnorm(q)
        #}

        nleqslv::nleqslv(2 * stats::qnorm(p), obj)$x
      }else{
        numeric(0)
      }
    }

    q0 <- rep(0, length(p[p == 0.5]) )
    if (length(p[p != 0.5]) < 1){
      q <- numeric(0)
    } else{
      q <- as.numeric(apply(matrix(p[p != 0.5], ncol = 1), 1, qtf))
    }

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

  ## Abbreviation
  out$name <- "Canonical slash"

  out
}

### Slash ----------------------------------------------------------------------
SL <- function(){
  out <- list()

  ## Incomplete gamma function
  ig <- function(a, x) gamma(a) * stats::pgamma(x, a)

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
          nu * x^(nu-1) * stats::pnorm(x * q)
        }

        stats::integrate(integrand, 0, 1)$val
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

        nleqslv::nleqslv(stats::qnorm(p) / (0.5^(1/nu)), obj)$x
      }else{
        numeric(0)
      }
    }

    q0 <- rep(0, length(p[p == 0.5]) )
    if (length(p[p != 0.5]) < 1){
      q <- numeric(0)
    } else{
      q <- as.numeric(apply(matrix(c(p[p != 0.5], nu[p != 0.5]), ncol = 2),
                            1, qtf))
    }

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

  ## Name
  out$name <- "Slash"

  out
}
