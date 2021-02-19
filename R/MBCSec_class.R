dmbcsec <- function(x, param, P = NULL, df = 4,
                    copula = c("gaussian", "st", "cauchy", "dexponential"),
                    gen = "NO"){

  ### Reading the copula generating density
  copula <- match.arg(copula, c("gaussian", "st", "cauchy", "dexponential"))

  ### Setting dimensions
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  n <- dim(x)[1]
  d <- dim(x)[2]

  ### Default association matrix
  if (is.null(P)) P <- diag(d)

  ### Marginal parameters
  mu <- param$mu
  sigma <- param$sigma
  lambda <- param$lambda
  nu <- param$nu

  ### Copula density argument
  if (length(gen) == 1){
    w <- pBCS(x, mu, sigma, lambda, nu, gen[1])
  }else{
    w <- matrix(apply(matrix(gen, ncol = 1), 1,
           function(gen) pBCS(x, mu, sigma, lambda, nu, gen))[, 1], ncol = 2)
  }

  if (copula != "st"){
    q <- BCSgen(ell(copula)$gen)$q(w)
    den <- log(ell(copula)$Md(q, P)) +
      apply(matrix(
        apply(matrix(gen, ncol = 1), 1,
        function(gen) log(dBCS(x, mu, sigma, lambda, nu, gen)))[, 1], ncol = 2
        ), 1, sum) -
      apply(matrix(log(BCSgen(ell(copula)$gen)$d(q^2)), ncol = d), 1, sum)

      #mvtnorm::dmvnorm(q, sigma = P, log = TRUE) +
      #apply(matrix(log(dBCS(x, mu, sigma, lambda, nu, gen[1])), ncol = d), 1, sum) -
      #apply(matrix(log(stats::dnorm(q)), ncol = d), 1, sum)
  } else {
    q <- BCSgen(ell(copula)$gen)$q(w, df)
    den <- log(ell(copula)$Md(q, P, df)) +
      apply(matrix(
        apply(matrix(gen, ncol = 1), 1,
        function(gen) log(dBCS(x, mu, sigma, lambda, nu, gen)))[, 1], ncol = 2
        ), 1, sum) -
      apply(matrix(log(BCSgen(ell(copula)$gen)$d(q^2, df)), ncol = d), 1, sum)

    #q <- stats::qt(w, df = df)
    #den <- mvtnorm::dmvt(q, delta = rep.int(0, d), sigma = P,
    #                     df = df, log = TRUE) +
    #apply(matrix(log(dBCS(x, mu, sigma, lambda, nu, gen[1])), ncol = d), 1, sum) -
    #apply(matrix(log(stats::dt(q, df = df)), ncol = d), 1, sum)

  }

  pmax(exp(den), .Machine$double.eps)

}
