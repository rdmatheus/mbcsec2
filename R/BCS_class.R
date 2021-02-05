#' @name BCS
#'
#' @title The Box-Cox Symmetric Class of Distributions
#'
#' @description Probability mass function, cumulative distribution function,
#' quantile function, and random generation for the Box-Cox symmetric (BCS)
#' class of distributions.
#'
#' @param x,q vector of non-negative quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param rep number of replicates with size \code{n} to return.
#' @param mu vector of location parameter values.
#' @param sigma vector of dispersion parameter values.
#' @param lambda vector of skewness parameter values.
#' @param nu vector of right tail/kurtosis parameter values. Not all
#'     distributions are indexed by this parameter. In this case,
#'     \code{nu = NULL}.
#' @param gen character specification passed as argument to \code{\link{BCSgen}}. It
#'     specifies the generating distribution for the BCS class of distributions.
#'     A table with the current available generating distributions can be seen
#'     in details.
#'
#' @details This set of functions represents the probability function, the
#'     cumulative distribution function, quantile function, and a random number
#'     generator for the Box-Cox symmetric class of distributions proposed by
#'     Ferrari and Fumes (2017). This class of distributions includes the
#'     Box-Cox t (Rigby and Stasinopoulos, 2006), Box-Cox Cole-Green
#'     (or Box-Cox normal; Cole and Green, 1992), Box-Cox power exponential
#'     (Rigby and Stasinopoulos, 2004) distributions, and the class of the
#'     log-symmetric distributions (Vanegas and Paula, 2016) as special cases.
#'     The current available generating distributions can be seen below.
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
#' @return \code{dBCS} returns the density function, \code{pBCS}
#' gives the distribution function, \code{qBCS} gives the quantile function,
#' and \code{rBCS} generates random observations.
#'
#' @references
#'  Cole, T., & Green, P.J. (1992). Smoothing reference centile curves: the LMS
#'      method and penalized likelihood. Stat. Med, 11, 1305--1319.
#'
#'  Ferrari, S. L., & Fumes, G. (2017). Box-Cox symmetric distributions and
#'      applications to nutritional data. AStA Advances in Statistical Analysis,
#'      101, 321--344.
#'
#'  Rigby, R. A., & Stasinopoulos, D. M. (2004). Smooth centile curves for skew
#'      and kurtotic data modelled using the Box-Cox power exponential
#'      distribution. Statistics in medicine, 23, 3053--3076.
#'
#'  Rigby, R. A., & Stasinopoulos, D. M. (2006). Using the Box-Cox t
#'      distribution in GAMLSS to model skewness and kurtosis. Statistical
#'      Modelling, 6, 209-229.
#'
#'  Vanegas, L. H., & Paula, G. A. (2016). Log-symmetric distributions:
#'      statistical properties and parameter estimation. Brazilian Journal of
#'      Probability and Statistics, 30, 196--220.
#'
#' @author Rodrigo M. R. Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' \dontrun{
#' ## Box-Cox t
#' y <- rBCS(500, 8, 0.15, 1, 10, gen = "ST")
#'
#' par(mfrow = c(1, 2))
#' hist(y, prob = TRUE, main = "Box Cox t")
#' curve(dBCS(x, 8, 0.15, 1, 10, gen = "ST"), add = TRUE, col = 2, lwd = 2)
#' plot(ecdf(y), main = "", xlab = "y")
#' curve(pBCS(x, 8, 0.15, 1, 10, gen = "ST"), add = TRUE, col = 2, lwd = 2)
#'
#' ## Box-Cox Cole-Green
#' y <- rBCS(500, 8, 0.15, 1, 10, gen = "NO")
#'
#' hist(y, prob = TRUE, main = "Box Cox Cole-Green")
#' curve(dBCS(x, 8, 0.15, 1, 10, gen = "NO"), add = TRUE, col = 2, lwd = 2)
#' plot(ecdf(y), main = "", xlab = "y")
#' curve(pBCS(x, 8, 0.15, 1, 10, gen = "NO"), add = TRUE, col = 2, lwd = 2)
#'
#' ## Box-Cox Cauchy
#' y <- rBCS(500, 8, 0.15, -1.2, gen = "CA")
#'
#' hist(y, prob = TRUE, main = "Box-Cox Cauchy")
#' curve(dBCS(x, 8, 0.15, -1.2, gen = "CA"), add = TRUE, col = 2, lwd = 2)
#' plot(ecdf(y), main = "", xlab = "y")
#' curve(pBCS(x, 8, 0.15, -1.2, gen = "CA"), add = TRUE, col = 2, lwd = 2)
#'
#' ## Box-Cox power exponential
#' y <- rBCS(500, 8, 0.15, 1, 10, gen = "PE")
#'
#' hist(y, prob = TRUE, main = "Box Cox PE")
#' curve(dBCS(x, 8, 0.15, 1, 10, gen = "PE"), add = TRUE, col = 2, lwd = 2)
#' plot(ecdf(y), main = "")
#' curve(pBCS(x, 8, 0.15, 1, 10, gen = "PE"), add = TRUE, col = 2, lwd = 2)
#' }
NULL

### Extended Box-Cox transformation
h <- function(y, mu, sigma, lambda){

  if (is.vector(y))
    y <- matrix(y, nrow = length(y))

  n <- dim(y)[1]
  d <- dim(y)[2]

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  lambda <- matrix(lambda, ncol = d)

  z <- matrix(NaN, n, d)

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))
  lambda <- do.call(rbind, replicate(n/dim(lambda)[1], lambda, simplify = FALSE))

  id1 <- which(y > 0 & mu > 0 & sigma > 0 & lambda != 0, arr.ind = TRUE)
  id2 <- which(y > 0 & mu > 0 & sigma > 0 & lambda == 0, arr.ind = TRUE)

  z[id1] <- ((y[id1]/mu[id1])^lambda[id1] - 1) /
    (sigma[id1] * lambda[id1])
  z[id2] <- log(y[id2] / mu[id2]) / sigma[id2]

  if(d == 1L) as.vector(z) else z
}

# Density
#' @rdname BCS
#' @export
dBCS <- function(x, mu, sigma, lambda, nu = NULL, gen = "NO"){

  # Extended Box-Cox transformation
  z <- h(x, mu, sigma, lambda)

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- dim(x)[1]
  d <- dim(x)[2]

  if (is.vector(z))
    z <- matrix(z, nrow = length(z))

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  lambda <- matrix(lambda, ncol = d)

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))
  lambda <- do.call(rbind, replicate(n/dim(lambda)[1], lambda, simplify = FALSE))

  if (!is.null(nu)){
    nu <- matrix(nu, ncol = d)
    nu <- do.call(rbind, replicate(n/dim(nu)[1], nu, simplify = FALSE))
  }

  # Generating distribution
  r <- BCSgen(gen)$d
  R <- BCSgen(gen)$p

  pmf <- matrix(0, n, d)

  # NaN index
  NaNid <- which(mu <= 0 | sigma <= 0, arr.ind = TRUE)

  pmf[NaNid] <- NaN

  # Positive density index
  id1 <- which(x > 0 & mu > 0 & sigma > 0 & lambda != 0 & !is.nan(pmf),
               arr.ind = TRUE)
  id2 <- which(x > 0 & mu > 0 & sigma > 0 & lambda == 0 & !is.nan(pmf),
               arr.ind = TRUE)


  pmf[id1] <- (x[id1]^(lambda[id1] - 1) * r(z[id1]^2, nu[id1])) /
    (R(1/(sigma[id1] * abs(lambda[id1])), nu[id1]) * sigma[id1] * mu[id1]^lambda[id1])

  pmf[id2] <- r(z[id2]^2, nu[id2]) / (sigma[id2] * x[id2])


  if(d == 1L) as.vector(pmf) else pmf

}

# Distribution function
#' @rdname BCS
#' @export
pBCS <- function(q, mu, sigma, lambda, nu = NULL, gen = "NO"){

  # Extended Box-Cox transformation
  z <- h(q, mu, sigma, lambda)

  if (is.vector(q))
    q <- matrix(q, nrow = length(q))

  n <- dim(q)[1]
  d <- dim(q)[2]

  if (is.vector(z))
    z <- matrix(z, nrow = length(z))

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  lambda <- matrix(lambda, ncol = d)

  distf <- matrix(NaN, n, d)

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify=FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify=FALSE))
  lambda <- do.call(rbind, replicate(n/dim(lambda)[1], lambda, simplify=FALSE))

  if (!is.null(nu)){
    nu <- matrix(nu, ncol = d)
    nu <- do.call(rbind, replicate(n/dim(nu)[1], nu, simplify = FALSE))
  }

  # Generating distribution
  r <- BCSgen(gen)$d
  R <- BCSgen(gen)$p

  id1 <- which(mu > 0 & sigma > 0 &lambda <= 0, arr.ind = TRUE)
  id2 <- which(mu > 0 & sigma > 0 &lambda > 0, arr.ind = TRUE)

  distf <- matrix(NA, n, d)
  distf[id1] <- R(z[id1], nu[id1]) / R(1 / (sigma[id1] * abs(lambda[id1])), nu[id1])
  distf[id2] <- (R(z[id2], nu[id2]) -  R(- 1 / (sigma[id2] * abs(lambda[id2])), nu[id2])) /
    R(1 / (sigma[id2] * abs(lambda[id2])), nu[id2])

  distf[which(q < 0, arr.ind = TRUE)] <- 0


  if (d == 1L) as.vector(distf) else distf
}

# Quantile function
#' @rdname BCS
#' @export
qBCS <- function(p, mu, sigma, lambda, nu = NULL, gen = "NO"){

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n <- dim(p)[1]
  d <- dim(p)[2]

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)
  lambda <- matrix(lambda, ncol = d)

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify=FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify=FALSE))
  lambda <- do.call(rbind, replicate(n/dim(lambda)[1], lambda, simplify=FALSE))

  if (!is.null(nu)){
    nu <- matrix(nu, ncol = d)
    nu <- do.call(rbind, replicate(n/dim(nu)[1], nu, simplify = FALSE))
  }

  # Generating distribution
  qr <- BCSgen(gen)$q
  R <- BCSgen(gen)$p

  q <- zp <- matrix(NaN, n, d)

  # Z_alpha
  zp[id1] <- qr(p[id1] * R(1 / (sigma[id1] * abs(lambda[id1])), nu[id1]), nu[id1])
  zp[id2] <- qr(1 - (1 - p[id2]) * R(1 / (sigma[id2] * abs(lambda[id2])), nu[id2]), nu[id2])

  id1 <- which(mu > 0 & sigma > 0 & lambda != 0, arr.ind = TRUE)
  id2 <- which(mu > 0 & sigma > 0 & lambda == 0, arr.ind = TRUE)

  q[id1] <- mu[id1] * (1 + sigma[id1] * lambda[id1] * zp[id1])^(1 / lambda[id1])
  q[id2] <- mu[id2] * exp(sigma[id2] * zp[id2])

  if (d == 1L) as.vector(q) else q
}

# Random generation
#' @rdname BCS
#' @export
rBCS <- function(n, mu, sigma, lambda, nu = NULL, gen = "NO", rep = 1L){

  if (rep == 1L){
    u <- stats::runif(n)
  }else{
    u <- matrix(stats::runif(n * rep), ncol = n)
  }

  x <- qBCS(u, mu, sigma, lambda, nu, gen)

  return(x)
}
