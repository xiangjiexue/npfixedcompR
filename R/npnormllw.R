npnormllw = R6::R6Class("npnormllw",
                       inherit = npnorm,
                       public = list(
                         beta = 1,
                         h = 10^-3,
                         w = NULL,
                         type = "npnorm",
                         initialize = function(data, mu0, pi0, beta, order = -3){
                           v = bin(data, order = order)
                           self$data = v$v
                           self$w = v$w
                           self$len = length(self$data)
                           self$h = 10^order
                           if (!missing(mu0)) self$mu0fixed = mu0
                           if (!missing(pi0)) self$pi0fixed = pi0
                           if (!missing(beta)) self$beta = beta
                           private$precompute = dnpdiscnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta, h = self$h)
                         },
                         modified = function(mu0, pi0, beta){
                           if (!missing(mu0)) self$mu0fixed = mu0
                           if (!missing(pi0)) self$pi0fixed = pi0
                           if (!missing(beta)) self$beta = beta
                           private$precompute = dnpdiscnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta, h = self$h)
                         },
                         lossfunction = function(mu0, pi0){
                           -sum(log(dnpdiscnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta, h = self$h) + private$precompute) * self$w) + sum(self$w) * log(self$h)
                         },
                         gradientfunction = function(mu, mu0, pi0, order = c(1, 0, 0)){
                           flexden = dnpdiscnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta, h = self$h)
                           fullden = flexden + private$precompute
                           murep = self$data - rep(mu, rep(self$len, length(mu)))
                           ans = vector("list", 3)
                           names(ans) = c("d0", "d1", "d2")
                           if (order[1] == 1){
                             temp = ddiscnorm(murep, sd = self$beta, h = self$h)
                             ans$d0 = .colSums((flexden - temp * sum(pi0)) / fullden * self$w, m = self$len, n = length(mu))
                           }
                           if (any(order[2:3] == 1)){
                             temp2 = dnorm(self$h + murep, sd = self$beta)
                             temp = temp2 - dnorm(murep, sd = self$beta)
                           }
                           if (order[2] == 1){
                             ans$d1 = .colSums(temp / fullden * self$w, m = self$len, n = length(mu)) * sum(pi0)
                           }
                           if (order[3] == 1){
                             temp1 = temp2 * self$h + temp * murep
                             ans$d2 = .colSums(temp1 / fullden * self$w, m = self$len, n = length(mu)) / self$beta^2
                           }

                           ans
                         },
                         computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                           mu0new = c(mu0, newweights)
                           pi0new = c(pi0, rep(0, length(newweights)))
                           sp = ddiscnorm(self$data, mean = rep(mu0new, rep(self$len, length(mu0new))),
                                          sd = self$beta, h = self$h)
                           dim(sp) = c(self$len, length(mu0new))
                           fp = drop(sp %*% pi0new) + private$precompute
                           S = sp / fp
                           a = 2 - private$precompute / fp
                           nw = pnnls(S * sqrt(self$w), a * sqrt(self$w), sum = 1 - sum(self$pi0fixed))$x
                           r = self$checklossfunction(mu0new, pi0new, nw - pi0new, crossprod(S, self$w))
                           self$collapsemix(r$pt, r$pr, tol)
                         },
                         estpi0dS = function(){
                           private$S1 = ddiscnorm(self$data, sd = self$beta, h = self$h)
                         },
                         estpi0d = function(mu0, pi0){
                           ans = vector("list", 2)
                           names(ans) = c("d2", "d3")
                           S = private$S1 / dnpdiscnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta, h = self$h) - 1
                           ans$d2 = -sum(S * self$w); ans$d3 = sum(S^2 * self$w)
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
                             self$result = list(iter = 0,
                                        family = self$type,
                                        max.gradient = self$gradientfunction(0, 0, 1, order = c(1, 0, 0))$d0,
                                        mix = list(pt = 0, pr = 1),
                                        beta = self$beta,
                                        ll = self$lossfunction(mu0 = 0, pi0 = 1),
                                        convergence = 0)
                           }else{
                             self$estpi0dS()
                             private$solveestpi0(init = dnpdiscnorm(0, mu0 = self$result$mix$pt, pi0 = self$result$mix$pr, sd = self$beta, h = self$h) / ddiscnorm(0, sd = self$beta, h = self$h),
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
makeobject.npnormllw = function(v, mu0, pi0, beta, order = -3){
  npnormllw$new(v, mu0, pi0, beta, order)
}
