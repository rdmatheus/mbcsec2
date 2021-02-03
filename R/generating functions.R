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
