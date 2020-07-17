# This file contains general functions used in npnorm family
# functions for normal mixture

#' The density and the distribution function of non-parametric normal distribution
#'
#' \code{dnpnorm} gives the density, \code{pnpnorm} gives the distribution function,
#' \code{pnpnorm1} focus on the more accurate but slower distribution function of
#' a non-parametric normal distribution
#'
#' @title non-parametric normal distribution
#' @param x vector of observations, vector of quantiles
#' @param mu0 the vector of support points
#' @param sd standard deviation.
#' @param pi0 the vector of weights correponding to the support points
#' @param lower.tail logical; if TRUE, the lower probability is computed
#' @param log,log.p logical; if TRUE, the result will be given in log scale.
#' @rdname npnorm
#' @export
dnpnorm = function(x, mu0 = 0, pi0 = 1, sd = 1, log = FALSE){
  # This version explicitly allow subprobability measures.
  # Hence no check on pi0
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = dnorm(x, mean = rep(mu0, rep(length(x), length(mu0))), sd = sd)
  dim(temp) = c(length(x), length(mu0))
  temp = drop(temp %*% pi0)
  if (log) log(temp) else temp
}

#' @rdname npnorm
#' @export
pnpnorm = function(x, mu0 = 0, pi0 = 1, sd = 1, lower.tail = TRUE, log.p = FALSE){
  # This version explicitly allow subprobability measures.
  # Hence no check on pi0
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = pnorm(x, mean = rep(mu0, rep(length(x), length(mu0))), sd = sd, lower.tail = lower.tail)
  dim(temp) = c(length(x), length(mu0))
  temp = drop(temp %*% pi0)
  if (log.p) log(temp) else temp
}

#' @rdname npnorm
#' @export
pnpnorm1 = function(x, mu0 = 0, pi0 = 1, sd = 1, lower.tail = TRUE, log.p = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  j0 = pi0 == 0
  if (sum(!j0) > 0){
    temp = pnorm(x, mean = rep(mu0[!j0], rep(length(x), sum(!j0))), sd = sd, lower.tail = lower.tail, log.p = TRUE) +
      rep(log(pi0[!j0]), rep(length(x), sum(!j0)))
    dim(temp) = c(length(x), sum(!j0))
    maxcoef = apply(temp, 1, max)
    temp = log(rowSums(exp(temp - maxcoef))) + maxcoef
  }else{
    temp = rep(-Inf, length(x))
  }

  if (log.p) temp else exp(temp)
}

# functions of discrete normal mixture

logspace.add = function(lx, ly){
  j0 = lx == -Inf & ly == -Inf
  ans = pmax(lx, ly) + log1p(exp(-abs(lx - ly)))
  ans[j0] = -Inf
  ans
}

log1mexp = function(x){
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

logspace.sub = function(lx, ly){
  lx + log1mexp(lx - ly)
}

#' The density and the distribution function of (non-parametric) discrete normal distribution
#'
#' \code{ddiscnorm} gives the density, \code{pdiscnorm} gives the distribution function of
#' the discrete normal distribution. \code{dnpdiscnorm} gives the density, \code{pnpdiscnorm}
#' gives the distribution function of the non-parametric discrete normal distribution.
#'
#' The function \code{pnpdiscnorm} uses \code{pnpnorm1} to compute the distribution
#' function.
#'
#' @title (non-parametric) discrete normal distribution
#' @param x vector of observations, vector of quantiles
#' @param mu0,mean the vector of support points
#' @param pi0 the vector of weights correponding to the support points \code{mu0}
#' @param sd standard deviation.
#' @param h the discretisation parameter.
#' @param lower.tail logical; if TRUE, the lower probability is computed
#' @param log,log.p logical; if TRUE, the result will be given in log scale.
#' @rdname discnorm
#' @export
ddiscnorm = function(x, mean = 0, sd = 1, h = 1, log = FALSE){
  LLL = max(length(x), length(mean))
  xx = rep(x, length.out = LLL)
  meanx = rep(mean, length.out = LLL)
  temp = logspace.sub(pnorm(xx + h, meanx, sd, log.p = TRUE), pnorm(xx, meanx, sd, log.p = TRUE))
  if (log) temp else exp(temp)
}

#' @rdname discnorm
#' @export
dnpdiscnorm = function(x, mu0 = 0, pi0 = 0, sd = 1, h = 1, lower.tail = TRUE, log.p = FALSE){
  j0 = pi0 == 0
  if (sum(!j0) > 0){
    temp = logspace.sub(pnorm(x + h, rep(mu0[!j0], rep(length(x), sum(!j0))), sd = sd, lower.tail = lower.tail, log.p = TRUE),
                        pnorm(x, rep(mu0[!j0], rep(length(x), sum(!j0))), sd = sd, lower.tail = lower.tail, log.p = TRUE)) +
      rep(log(pi0[!j0]), rep(length(x), sum(!j0)))
    dim(temp) = c(length(x), sum(!j0))
    maxcoef = apply(temp, 1, max)
    temp = log(rowSums(exp(temp - maxcoef))) + maxcoef
  }else{
    temp = rep(-Inf, length(x))
  }

  if (log.p) temp else exp(temp)
}

#' @rdname discnorm
#' @export
pdiscnorm = function(x, mean =  0, sd = 1, h = 1, lower.tail = TRUE, log.p = FALSE){
  pnorm(x + h, mean = mean, sd = sd, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname discnorm
#' @export
pnpdiscnorm = function(x, mu0 = 0, pi0 = 0, sd = 1, h = 1, lower.tail = TRUE, log.p = FALSE){
  pnpnorm1(x + h, mu0 = mu0, pi0 = pi0, sd = sd, lower.tail = lower.tail, log.p = log.p)
}

npnorm = R6::R6Class("npnorm",
                 inherit = npfixedcompR,
                 public = list(
                   setgridpoints = function(grid=100) {
                     rx = range(self$data)
                     breaks = pmax(ceiling(diff(rx) / (5*self$beta)), 5)   # number of breaks
                     if (is.null(self$w)) {w = 1} else {w = self$w}
                     r = whist(self$data, w, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
                     i = r$density != 0
                     i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
                     m = sum(i)
                     k = pmax(ceiling(grid / m), 10)           # at least 10 in each interval
                     d = r$breaks[2] - r$breaks[1]
                     s = r$breaks[-length(r$breaks)][i]
                     private$gridpoints = sort(c(rx[1], rep(s, rep(k,m)) + d * (1:k-0.5)/k, rx[2]), decreasing = FALSE)
                   }
                 ))



