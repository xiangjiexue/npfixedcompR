npnormll = R6::R6Class("npnormll",
                   inherit = npnorm,
                   public = list(
                     beta = 1,
                     type = "npnorm",
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
                       if (any(self$compareattr(mu0, pi0))){
                         self$setflexden(mu0, pi0)
                       }
                       -sum(log(private$flexden$fullden))
                     },
                     setflexden = function(mu0, pi0){
                       r = self$compareattr(mu0, pi0)
                       if (r[1]){
                         private$flexden$dens = dnorm(self$data, mean = rep(mu0, rep(self$len, length(mu0))), sd = self$beta)
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
                       if (any(self$compareattr(mu0, pi0))){
                         self$setflexden(mu0, pi0)
                       }
                       mu0new = c(mu0, newweights)
                       pi0new = c(pi0, rep(0, length(newweights)))
                       sp = cbind(private$flexden$dens, matrix(dnorm(self$data, mean = rep(newweights, rep(self$len, length(newweights))), sd = self$beta), nrow = self$len, ncol = length(newweights)))
                       fp = private$flexden$fullden
                       S = sp / fp
                       a = 2 - private$precompute / fp
                       nw = pnnls(S, a, sum = 1 - sum(self$pi0fixed))$x
                       r = self$checklossfunction(mu0new, pi0new, nw - pi0new, colSums(S))
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
                                  family = x$type,
                                  max.gradient = self$gradientfunction(0, 0, 1, order = c(1, 0, 0))$d0,
                                  mix = list(pt = 0, pr = 1),
                                  beta = self$beta,
                                  ll = self$lossfunction(mu0 = 0, pi0 = 1),
                                  convergence = 0)
                       }else{
                         self$estpi0dS()
                         private$solveestpi0(init = dnpnorm(0, mu0 = self$result$mix$pt, pi0 = self$result$mix$pr, sd = self$beta) * sqrt(2 * base::pi) * self$beta,
                                         val = -r0ll - val, tol = tol, maxiter = maxiter, verbose = verbose)
                       }
                     }
                   ),
                   private = list(
                     precompute = NULL,
                     flexden = NULL,
                     mflag = "d1",
                     S1 = NULL
                   ))

#' @rdname makeobject
#' @export
makeobject.npnormll = function(v, mu0, pi0, beta){
  npnormll$new(v, mu0, pi0, beta)
}
