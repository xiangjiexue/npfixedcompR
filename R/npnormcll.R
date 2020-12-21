dnormc = function(x, mean = 0, n, log = FALSE){
  if (any(abs(x) > 1) | any(abs(mean) > 1))
    stop("Error in specifying data or mean")
  dnorm(x, mean = mean, sd = (1 - mean^2) / sqrt(n), log = log)
}

pnormc = function(x, mean = 0, n, lower.tail = TRUE, log.p = FALSE){
  if (any(abs(x) > 1) | any(abs(mean) > 1))
    stop("Error in specifying data or mean")
  pnorm(x, mean = mean, sd = (1 - mean^2) / sqrt(n), lower.tail = lower.tail, log.p = log.p)
}

dnpnormc = function(x, mu0 = 0, pi0 = 0, n, log = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = dnormc(x, mean = rep(mu0, rep(length(x), length(mu0))), n = n)
  dim(temp) = c(length(x), length(mu0))
  temp = drop(temp %*% pi0)
  if (log) log(temp) else temp
}

pnpnormc = function(x, mu0 = 0, pi0 = 0, n, lower.tail = TRUE, log.p = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = pnormc(x, mean = rep(mu0, rep(length(x), length(mu0))), n = n, lower.tail = lower.tail)
  dim(temp) = c(length(x), length(mu0))
  temp = drop(temp %*% pi0)
  if (log.p) log(temp) else temp
}

npnormcll = R6::R6Class("npnormcll",
                    inherit = npfixedcompR,
                    public = list(
                      setgridpoints = function(grid=100){
                        private$gridpoints = seq(from = min(self$data), to = max(self$data), length.out = grid)
                      },
                      initpoints = function(){
                        rx = range(self$data)
                        breaks = pmax(ceiling(diff(rx) / (5 / sqrt(self$beta))), 10)   # number of breaks
                        r = whist(self$data, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
                        r$density = pmax(0, r$density  / sum(r$density) - pnpnormc(r$breaks[-1], mu0 = self$mu0fixed, pi0 = self$pi0fixed, n = self$beta) +
                                           pnpnormc(r$breaks[-length(r$breaks)], mu0 = self$mu0fixed, pi0 = self$pi0fixed, n = self$beta))
                        list(pt = r$mids[r$density != 0], pr = r$density[r$density != 0] / sum(r$density))
                      },
                      beta = NULL,
                      type = "npnormc",
                      initialize = function(data, mu0, pi0, beta){
                        self$data = data
                        self$len = length(data)
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$precompute = dnpnormc(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, n = self$beta)
                      },
                      modified = function(mu0, pi0, beta){
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$precompute = dnpnormc(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, n = self$beta)
                      },
                      lossfunction = function(mu0, pi0){
                        -sum(log(dnpnormc(self$data, mu0 = mu0, pi0 = pi0, n = self$beta) + private$precompute))
                      },
                      gradientfunction = function(mu, mu0, pi0, order = c(1, 0, 0)){
                        flexden = dnpnormc(self$data, mu0 = mu0, pi0 = pi0, n = self$beta)
                        fullden = flexden + private$precompute
                        ans = vector("list", 3)
                        names(ans) = c("d0", "d1", "d2")
                        if (length(mu) > 0){
                          murep = rep(mu, rep(self$len, length(mu)))
                          temp = dnormc(self$data, mean = murep, n = self$beta)
                          # only d0 is implemented
                          if (order[1] == 1){
                            ans$d0 = .colSums((flexden - temp * sum(pi0)) / fullden, m = self$len, n = length(mu))
                          }
                          if (order[2] == 1){
                            temp2 = ((self$beta + 4) * murep^3 - 3 * self$beta * murep^2 * self$data + murep * (2 * self$beta * self$data^2 + self$beta - 2) -
                                       self$beta * self$data - 2 * murep^5) / (1 - murep^2)^3
                            ans$d1 = .colSums(temp2 * temp / fullden, m = self$len, n = length(mu))  * sum(pi0)
                          }

                        }

                        ans
                      },
                      computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                        mu0new = c(mu0, newweights)
                        pi0new = c(pi0, rep(0, length(newweights)))
                        sp = dnormc(self$data, mean = rep(mu0new, rep(self$len, length(mu0new))),
                                   n = self$beta)
                        dim(sp) = c(self$len, length(mu0new))
                        fp = drop(sp %*% pi0new) + private$precompute
                        S = sp / fp
                        a = 2 - private$precompute / fp
                        nw = pnnls(S, a, sum = 1 - sum(self$pi0fixed))$x
                        r = self$checklossfunction(mu0new, pi0new, nw - pi0new, colSums(S))
                        self$collapsemix(r$pt, r$pr, tol)
                      },
                      estpi0 = function(val = -log(0.5), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
                        stop("It is not possible to estimate the proportion of null.")
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
makeobject.npnormcll = function(v, mu0, pi0, beta, order){
  npnormcll$new(v, mu0, pi0, beta)
}
