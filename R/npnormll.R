npnormll = R6::R6Class("npnormll",
                   inherit = npnorm,
                   public = list(
                     beta = 1,
                     type = "npnormll",
                     initialize = function(data, mu0, pi0, beta){
                        self$data = data
                        self$len = length(data)
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$flexden = list(flexden = numeric(self$len), fullden = numeric(self$len))
                        attr(private$flexden, "mu0") = NULL
                        attr(private$flexden, "pi0") = NULL
                        private$precompute = dnpnorm(data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta)
                     },
                     modified = function(mu0, pi0, beta){
                       if (!missing(mu0)) self$mu0fixed = mu0
                       if (!missing(pi0)) self$pi0fixed = pi0
                       if (!missing(beta)) self$beta = beta
                       private$precompute = dnpnorm(data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta)
                     },
                     lossfunction = function(mu0, pi0){
                       -sum(log(dnpnorm(self$data, mu0 = mu0, pi0 = pi0, self$beta) + private$precompute))
                     },
                     setflexden = function(mu0, pi0){
                       temp = dnpnorm(self$data, mu0, pi0, self$beta)
                       private$flexden$flexden = temp
                       private$flexden$fullden = temp + private$precompute
                       attr(private$flexden, "mu0") = mu0
                       attr(private$flexden, "pi0") = pi0
                     },
                     gradientfunction = function(mu, mu0, pi0, order = c(1, 0, 0)){
                       if (self$compareattr(mu0, pi0)){
                         self$setflexden(mu0, pi0)
                       }
                       murep = self$data - rep(mu, rep(self$len, length(mu)))
                       temp = dnorm(murep, sd = self$beta) * sum(pi0)
                       ans = vector("list", 3)
                       names(ans) = c("d0", "d1", "d2")
                       if (order[1] == 1){
                         ans$d0 = .colSums((private$flexden$flexden - temp) / private$flexden$fullden,
                                           m = self$len, n = length(mu))
                       }
                       if (any(order[2:3] == 1)){
                         temp2 = temp / private$flexden$fullden
                       }
                       if (order[2] == 1){
                         ans$d1 = .colSums(temp2 * murep, m = self$len, n = length(mu)) / -self$beta^2
                       }
                       if (order[3] == 1){
                         ans$d2 = .colSums(temp2 * (murep^2 - self$beta^2), m = self$len, n = length(mu)) / -self$beta^4
                       }

                       ans
                     },
                     computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                       mu0new = c(mu0, newweights)
                       pi0new = c(pi0, rep(0, length(newweights)))
                       sp = dnorm(self$data, mean = rep(mu0new, rep(self$len, length(mu0new))), sd = self$beta)
                       fp = .rowSums(sp[1:(self$len * length(mu0))] * rep(pi0, rep(self$len, length(pi0))), m = self$len, n = length(pi0)) + private$precompute
                       S = sp / fp
                       dim(S) = c(self$len, length(pi0new))
                       a = 2 - private$precompute / fp
                       nw = pnnls(S, a, sum = 1 - sum(self$pi0fixed))$x
                       r = self$checklossfunction(mu0new, pi0new, nw - pi0new, colSums(S), tol)
                       self$collapsemix(r$pt, r$pr, tol)
                     },
                     estpi0dS = function(){
                       private$S1 = dnorm(self$data, sd = self$beta)
                     },
                     estpi0d = function(mu0, pi0){
                       ans = vector("list", 2)
                       names(ans) = c("d2", "d3")
                       S = private$S1 / dnpnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta) - 1
                       ans$d2 = -sum(S); ans$d3 = sum(S^2)
                       ans
                     }
                   ),
                   private = list(
                     precompute = NULL,
                     flexden = NULL,
                     methodflag = "d2",
                     S1 = NULL
                   ))

#' @rdname makeobject
#' @export
makeobject.npnormll = function(v, mu0, pi0, beta){
  npnormll$new(v, mu0, pi0, beta)
}

#' @rdname estpi0
#' @export
estpi0.npnormll = function(x, val = 0.5 * log(x$len), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x$modified(pi0 = 1 - tol / 2)
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x$modified(pi0 = 0)
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll - r0$ll < val){
    r = list(iter = 0,
             family = "npnorm",
             max.gradient = x$gradientfunction(0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             beta = x$beta,
             ll = x$lossfunction(mu0 = 0, pi0 = 1),
             convergence = 0)
  }else{
    x$estpi0dS()
    r = solveestpi0(x = x, init = dnpnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$beta) * sqrt(2 * base::pi) * x$beta,
                    val = -r0$ll - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
