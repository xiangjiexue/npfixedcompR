
#' The density and the distribution function of non-parametric t distribution
#'
#' \code{dnpnorm} gives the density, \code{pnpnorm} gives the distribution function,
#'
#' @title non-parametric t distribution
#' @param x vector of observations, vector of quantiles
#' @param mu0 the vector of support points
#' @param df degree of freedom.
#' @param pi0 the vector of weights correponding to the support points
#' @param lower.tail logical; if TRUE, the lower probability is computed
#' @param log,log.p logical; if TRUE, the result will be given in log scale.
#' @rdname npt
#' @export
dnpt = function(x, mu0 = 0, pi0 = 0, df, log = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = .rowSums(dt(x, ncp = rep(mu0, rep(length(x), length(mu0))), df = df) *
                    rep(pi0, rep(length(x), length(pi0))), m = length(x), n = length(mu0))
  if (log) log(temp) else temp
}

#' @rdname npt
#' @export
pnpt = function(x, mu0 = 0, pi0 = 0, df, lower.tail = TRUE, log.p = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = .rowSums(pt(x, ncp = rep(mu0, rep(length(x), length(mu0))), df = df, lower.tail = lower.tail) *
                    rep(pi0, rep(length(x), length(pi0))), m = length(x), n = length(mu0))
  if (log.p) log(temp) else temp
}

nptll = R6::R6Class("nptll",
                    inherit = npfixedcompR,
                    public = list(
                      setgridpoints = function(grid=100){
                        rx = range(self$data)
                        if (is.finite(self$beta) & self$beta >= 3) {
                          fac = sqrt(self$beta / (self$beta - 2))
                          if (rx[2] > 0) {
                            rx[2] = rx[2] * sqrt((2 * self$beta + 5) / 2 / self$beta)
                          } else {
                            rx[2] = rx[2] * sqrt((self$beta + 1) / self$beta)
                          }

                          if (rx[1] > 0) {
                            rx[1] = rx[1] * sqrt((self$beta + 1) / self$beta)
                          } else {
                            rx[1] = rx[1] * sqrt((2 * self$beta + 5) / 2 / self$beta)
                          }
                        } else {fac = 1}
                        breaks = pmax(ceiling(diff(rx) / (5 * fac)), 5)   # number of breaks
                        r = hist(self$data, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
                        i = r$density != 0
                        i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
                        m = sum(i)
                        k = pmax(ceiling(grid / m), 10)           # at least 10 in each interval
                        d = r$breaks[2] - r$breaks[1]
                        s = r$breaks[-length(r$breaks)][i]
                        private$gridpoints = sort(c(rx[1], rep(s, rep(k,m)) + d * (1:k-0.5)/k, rx[2]), decreasing = FALSE)
                      },
                      initpoints = function(){
                        rx = range(self$data)
                        if (is.finite(self$beta) & self$beta >= 3) {fac = sqrt(self$beta / (self$beta - 2))} else fac = 1
                        breaks = pmax(ceiling(diff(rx) / (5 * fac)), 10)   # number of breaks
                        r = hist(self$data, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
                        r$density = pmax(0, r$density - pnpt(r$breaks[-1], mu0 = self$mu0fixed, pi0 = self$pi0fixed, df = self$beta) +
                                           pnpt(r$breaks[-length(r$breaks)], mu0 = self$mu0fixed, pi0 = self$pi0fixed, df = self$beta))
                        list(pt = r$mids[r$density != 0], pr = r$density[r$density != 0] / sum(r$density))
                      },
                      beta = Inf,
                      type = "npt",
                      initialize = function(data, mu0, pi0, beta){
                        self$data = data
                        self$len = length(data)
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$flexden = list(flexden = numeric(self$len), fullden = numeric(self$len))
                        attr(private$flexden, "mu0") = NULL
                        attr(private$flexden, "pi0") = NULL
                        private$precompute = dnpt(data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, df = self$beta)
                      },
                      modified = function(mu0, pi0, beta){
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$precompute = dnpt(data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, df = self$beta)
                      },
                      lossfunction = function(mu0, pi0){
                        if (any(self$compareattr(mu0, pi0))){
                          self$setflexden(mu0, pi0)
                        }
                        -sum(log(private$flexden$fullden))
                      },
                      setflexden = function(mu0, pi0){
                        r = self$compareattr(mu0, pi0)
                        if (r[1]){
                          private$flexden$dens = dt(self$data, ncp = rep(mu0, rep(self$len, length(mu0))), df = self$beta)
                          dim(private$flexden$dens) = c(self$len, length(mu0))
                        }
                        temp = drop(private$flexden$dens %*% pi0)
                        private$flexden$flexden = temp
                        private$flexden$fullden = temp + private$precompute
                        attr(private$flexden, "mu0") = mu0
                        attr(private$flexden, "pi0") = pi0
                      },
                      gradientfunction = function(mu, mu0, pi0, order = c(1, 0, 0)){
                        if (any(self$compareattr(mu0, pi0))){
                          self$setflexden(mu0, pi0)
                        }
                        temp = dt(self$data, ncp = rep(mu, rep(self$len, length(mu))), df = self$beta) * sum(pi0)
                        ans = vector("list", 3)
                        names(ans) = c("d0", "d1", "d2")

                        # only d0 is implemented
                        if (order[1] == 1){
                          ans$d0 = .colSums((private$flexden$flexden - temp) / private$flexden$fullden,
                                            m = self$len, n = length(mu))
                        }

                        ans
                      },
                      computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                        if (any(self$compareattr(mu0, pi0))){
                          self$setflexden(mu0, pi0)
                        }
                        mu0new = c(mu0, newweights)
                        pi0new = c(pi0, rep(0, length(newweights)))
                        sp = cbind(private$flexden$dens, matrix(dt(self$data, ncp = rep(newweights, rep(self$len, length(newweights))), df = self$beta), nrow = self$len, ncol = length(newweights)))
                        fp = private$flexden$fullden
                        S = sp / fp
                        a = 2 - private$precompute / fp
                        nw = pnnls(S, a, sum = 1 - sum(self$pi0fixed))$x
                        r = self$checklossfunction(mu0new, pi0new, nw - pi0new, colSums(S))
                        self$collapsemix(r$pt, r$pr, tol)
                      },
                      estpi0dS = function(){
                        private$S1 = dt(self$data, df = self$beta)
                      },
                      estpi0d = function(mu0, pi0){
                        ans = vector("list", 2)
                        names(ans) = c("d2", "d3")
                        S = private$S1 / dnpt(self$data, mu0 = mu0, pi0 = pi0, df = self$beta) - 1
                        ans$d2 = -sum(S); ans$d3 = sum(S^2)
                        ans
                      },
                      estpi0 = function(val = -log(0.5), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
                        self$modified(pi0 = 1 - tol / 2)
                        self$computemixdist(mix = mix, tol = tol, maxiter = maxiter)
                        r1ll = self$result$ll
                        self$modified(pi0 = 0)
                        self$computemixdist(mix = mix, tol = tol, maxiter = maxiter)
                        r0ll = self$result$ll

                        if (r1ll - r0ll < val){
                          r = list(iter = 0,
                                   family = self$type,
                                   max.gradient = self$gradientfunction(0, 0, 1, order = c(1, 0, 0))$d0,
                                   mix = list(pt = 0, pr = 1),
                                   beta = self$beta,
                                   ll = self$lossfunction(mu0 = 0, pi0 = 1),
                                   convergence = 0)
                        }else{
                          self$estpi0dS()
                          private$solveestpi0(init = dnpt(0, mu0 = self$result$mix$pt, pi0 = self$result$mix$pr, df = self$beta) / dt(0, df = self$beta),
                                           val = -r0ll - val, tol = tol, maxiter = maxiter, verbose = verbose)
                        }
                      }
                    ),
                    private = list(
                      precompute = NULL,
                      flexden = NULL,
                      mflag = "d0",
                      S1 = NULL
                    ))

#' @rdname makeobject
#' @export
makeobject.nptll = function(v, mu0, pi0, beta){
  nptll$new(v, mu0, pi0, beta)
}
