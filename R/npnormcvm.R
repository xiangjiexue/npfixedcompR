npnormcvm = R6::R6Class("npnormcvm",
                       inherit = npnorm,
                       public = list(
                         beta = 1,
                         type = "npnorm",
                         initialize = function(data, mu0, pi0, beta){
                           self$data = sort(data, decreasing = FALSE)
                           self$len = length(data)
                           if (!missing(mu0)) self$mu0fixed = mu0
                           if (!missing(pi0)) self$pi0fixed = pi0
                           if (!missing(beta)) self$beta = beta
                           self$estpi0dS()
                           private$flexden = list(flexden = numeric(self$len), fullden = numeric(self$len))
                           attr(private$flexden, "mu0") = NULL
                           attr(private$flexden, "pi0") = NULL
                           private$precompute = private$S1$a2 - pnpnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta)
                         },
                         modified = function(mu0, pi0, beta){
                           if (!missing(mu0)) self$mu0fixed = mu0
                           if (!missing(pi0)) self$pi0fixed = pi0
                           if (!missing(beta)) self$beta = beta
                           private$precompute = private$S1$a2 - pnpnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta)
                         },
                         lossfunction = function(mu0, pi0){
                           sum((pnpnorm(self$data, mu0, pi0, self$beta) - private$precompute)^2)
                         },
                         setflexden = function(mu0, pi0){
                           temp = pnpnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta)
                           private$flexden$flexden = temp
                           private$flexden$fullden = temp - private$precompute
                           attr(private$flexden, "mu0") = mu0
                           attr(private$flexden, "pi0") = pi0
                         },
                         gradientfunction = function(mu, mu0, pi0, order = c(1, 0, 0)){
                           if (self$compareattr(mu0, pi0)){
                             self$setflexden(mu0, pi0)
                           }
                           murep = self$data - rep(mu, rep(self$len, length(mu)))
                           ans = vector("list", 3)
                           names(ans) = c("d0", "d1", "d2")
                           if (order[1] == 1){
                             temp = pnorm(murep, sd = self$beta) * sum(pi0)
                             ans$d0 = .colSums((temp - private$flexden$flexden) * private$flexden$fullden, m = self$len, n = length(mu))
                           }
                           if (any(order[2:3] == 1)){
                             temp = dnorm(murep, sd = self$beta) * private$flexden$fullden
                           }
                           if (order[2] == 1){
                             ans$d1 = .colSums(temp, m = self$len, n = length(mu)) * -2 * sum(pi0)
                           }
                           if (order[3] == 1){
                             ans$d2 = .colSums(temp * murep, m = self$len, n = length(mu)) * -2 * sum(pi0)
                           }

                           ans
                         },
                         computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                           mu0new = c(mu0, newweights)
                           S = pnorm(self$data, mean = rep(mu0new, rep(self$len, length(mu0new))), sd = self$beta)
                           dim(S) = c(self$len, length(mu0new))
                           pi0new = pnnls(S, private$precompute, sum = 1 - sum(x$pi0fixed))$x
                           self$collapsemix(mu0new, pi0new, tol)
                         },
                         estpi0dS = function(){
                           private$S1 = list(a1 = pnorm(self$data, sd = self$beta),
                                             a2 = (2 * (1 : self$len) - 1) / 2 / self$len)
                         },
                         estpi0d = function(mu0, pi0){
                           ans = vector("list", 2)
                           names(ans) = c("d2", "d3")
                           pn = pnpnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta)
                           S = private$S1$a1 - pn
                           ans$d2 = 2 * sum(S * (pn - private$S1$a2)); ans$d3 = 2 * sum(S^2)
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
makeobject.npnormcvm = function(v, mu0, pi0, beta){
  npnormcvm$new(v, mu0, pi0, beta)
}

#' @rdname estpi0
#' @export
estpi0.npnormcvm = function(x, val = qCvM(0.5, lower.tail = FALSE), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x$modified(pi0 = 1 - tol / 2)
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x$modified(pi0 = 0)
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll + 1/ 12 / x$len < val){
    r = list(iter = 0,
             family = x$type,
             max.gradient = x$gradientfunction(0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             ll = x$lossfunction(mu0 = 0, pi0 = 1),
             beta = x$beta,
             convergence = 0)
  }else{
    x$estpi0dS()
    r = solveestpi0(x = x, init = dnpnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$beta) * sqrt(2 * base::pi) * x$beta,
                    val = 1/ 12 / x$len - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
