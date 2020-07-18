npnormad = R6::R6Class("npnormad",
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
                           private$flexden = list(flexden = numeric(self$len), flexden2 = numeric(self$len),
                                                  fullden = numeric(self$len), fullden2 = numeric(self$len),
                                                  addconst = numeric(1))
                           attr(private$flexden, "mu0") = NULL
                           attr(private$flexden, "pi0") = NULL
                           a1 = seq(from = 1 / self$len, to = (2 * self$len - 1) / self$len, length = self$len)
                           private$precompute = list(a1 = a1, a2 = rev(a1),
                                                     precompute1 = pnpnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta),
                                                     precompute2 = pnpnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta, lower.tail = FALSE))
                         },
                         modified = function(mu0, pi0, beta){
                           if (!missing(mu0)) self$mu0fixed = mu0
                           if (!missing(pi0)) self$pi0fixed = pi0
                           if (!missing(beta)) self$beta = beta
                           private$precompute$precompute1 = pnpnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta)
                           private$precompute$precompute2 = pnpnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta, lower.tail = FALSE)
                         },
                         lossfunction = function(mu0, pi0){
                           if (any(self$compareattr(mu0, pi0))){
                             self$setflexden(mu0, pi0)
                           }
                           -sum(log(private$flexden$fullden) * private$precompute$a1 +
                                  log(private$flexden$fullden2) * private$precompute$a2)
                         },
                         setflexden = function(mu0, pi0){
                           r = self$compareattr(mu0, pi0)
                           if (r[1]){
                             private$flexden$dens = pnorm(self$data, mean = rep(mu0, rep(self$len, length(mu0))), sd = self$beta)
                             dim(private$flexden$dens) = c(self$len, length(mu0))
                             private$flexden$dens2 = pnorm(self$data, mean = rep(mu0, rep(self$len, length(mu0))), sd = self$beta, lower.tail = FALSE)
                             dim(private$flexden$dens2) = c(self$len, length(mu0))
                           }
                           temp = drop(private$flexden$dens %*% pi0)
                           temp2 = drop(private$flexden$dens2 %*% pi0)
                           private$flexden$flexden = temp
                           private$flexden$fullden = temp + private$precompute$precompute1
                           private$flexden$flexden2 = temp2
                           private$flexden$fullden2 = temp2 + private$precompute$precompute2
                           private$flexden$addconst = -sum(private$precompute$precompute1 * private$precompute$a1 / private$flexden$fullden +
                                                             private$precompute$precompute2 * private$precompute$a2 / private$flexden$fullden2) +
                             2 * self$len
                           attr(private$flexden, "mu0") = mu0
                           attr(private$flexden, "pi0") = pi0
                         },
                         gradientfunction = function(mu, mu0, pi0, order = c(1, 0, 0)){
                           if (any(self$compareattr(mu0, pi0))){
                             self$setflexden(mu0, pi0)
                           }
                           murep = self$data - rep(mu, rep(self$len, length(mu)))
                           ans = vector("list", 3)
                           names(ans) = c("d0", "d1", "d2")
                           if (order[1] == 1){
                             ans$d0 = .colSums(pnorm(murep, sd = self$beta) * private$precompute$a1 / private$flexden$fullden +
                                                 pnorm(murep, sd = self$beta, lower.tail = FALSE) * private$precompute$a2 / private$flexden$fullden2,
                                               m = self$len, n = length(mu)) * -sum(pi0) + private$flexden$addconst
                           }
                           if (any(order[2:3] == 1)){
                             temp = dnorm(murep, sd = self$beta) * (private$precompute$a1 / private$flexden$fullden -
                                                                      private$precompute$a2 / private$flexden$fullden2)
                           }
                           if (order[2] == 1){
                             ans$d1 = .colSums(temp, m = self$len, n = length(mu)) * sum(pi0)
                           }
                           if (order[3] == 1){
                             ans$d2 = .colSums(murep * temp, m = self$len, n = length(mu)) * sum(pi0) / self$beta^2
                           }

                           ans
                         },
                         computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                           if (any(self$compareattr(mu0, pi0))){
                             self$setflexden(mu0, pi0)
                           }
                           mu0new = c(mu0, newweights)
                           pi0new = c(pi0, rep(0, length(newweights)))
                           sf = cbind(private$flexden$dens, matrix(pnorm(self$data, mean = rep(newweights, rep(self$len, length(newweights))), sd = self$beta), nrow = self$len, ncol = length(newweights)))
                           ss = private$flexden$fullden
                           S = sf / ss

                           uf = cbind(private$flexden$dens2, matrix(pnorm(self$data, mean = rep(newweights, rep(self$len, length(newweights))), sd = self$beta, lower.tail = FALSE), nrow = self$len, ncol = length(newweights)))
                           us = private$flexden$fullden2
                           U = uf / us

                           Sma1 = S * private$precompute$a1
                           Uma2 = U * private$precompute$a2
                           S1 = crossprod(S, Sma1) + crossprod(U, Uma2)
                           S2 = colSums(Sma1) + colSums(Uma2)

                           nw = pnnqp(S1, -2 * S2 + crossprod(S, private$precompute$precompute1 * private$precompute$a1 / ss) +
                                               crossprod(U, private$precompute$precompute2 * private$precompute$a2 / us), sum = 1 - sum(self$pi0fixed))$x
                           r = self$checklossfunction(mu0new, pi0new, nw - pi0new, S2)
                           self$collapsemix(r$pt, r$pr, tol)
                         },
                         estpi0dS = function(){
                           private$S1 = list(S1 = pnorm(self$data, sd = self$beta),
                                             S2 = pnorm(self$data, sd = self$beta, lower = FALSE))
                         },
                         estpi0d = function(mu0, pi0){
                           ans = vector("list", 2)
                           names(ans) = c("d2", "d3")
                           S1 = private$S1$S1 / pnpnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta) - 1
                           S2 = private$S1$S2 / pnpnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta, lower.tail = FALSE) - 1
                           ans$d2 = -sum(S1 * private$precompute$a1 + S2 * private$precompute$a2)
                           ans$d3 = sum(S1^2 * private$precompute$a1 + S2^2 * private$precompute$a2)
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
makeobject.npnormad = function(v, mu0, pi0, beta){
  npnormad$new(v, mu0, pi0, beta)
}

#' @rdname estpi0
#' @export
estpi0.npnormad = function(x, val = qAD(0.5, lower.tail = FALSE), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x$modified(pi0 = 1 - tol / 2)
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x$modified(pi0 = 0)
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll - x$len < val){
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
                    val = -x$len - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
