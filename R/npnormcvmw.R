npnormcvmw = R6::R6Class("npnormcvmw",
                        inherit = npnorm,
                        public = list(
                          beta = 1,
                          type = "npnorm",
                          h = 10^-3,
                          w = NULL,
                          initialize = function(data, mu0, pi0, beta, order = -3){
                            v = bin(data, order = order)
                            self$data = v$v
                            self$w = v$w
                            self$len = length(self$data)
                            self$h = 10^order
                            if (!missing(mu0)) self$mu0fixed = mu0
                            if (!missing(pi0)) self$pi0fixed = pi0
                            if (!missing(beta)) self$beta = beta
                            self$estpi0dS()
                            private$flexden = list(flexden = numeric(self$len), fullden = numeric(self$len))
                            attr(private$flexden, "mu0") = NULL
                            attr(private$flexden, "pi0") = NULL
                            private$precompute = private$S1$a2 - pnpdiscnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta, h = self$h)
                          },
                          modified = function(mu0, pi0, beta){
                            if (!missing(mu0)) self$mu0fixed = mu0
                            if (!missing(pi0)) self$pi0fixed = pi0
                            if (!missing(beta)) self$beta = beta
                            private$precompute = private$S1$a2 - pnpdiscnorm(self$data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, sd = self$beta, h = self$h)
                          },
                          lossfunction = function(mu0, pi0){
                            if (any(self$compareattr(mu0, pi0))){
                              self$setflexden(mu0, pi0)
                            }
                            sum(private$flexden$fullden^2 * self$w)
                          },
                          setflexden = function(mu0, pi0){
                            r = self$compareattr(mu0, pi0)
                            if (r[1]){
                              private$flexden$dens = pdiscnorm(self$data, mean = rep(mu0, rep(self$len, length(mu0))), sd = self$beta, h = self$h)
                              dim(private$flexden$dens) = c(self$len, length(mu0))
                            }
                            temp = drop(private$flexden$dens %*% pi0)
                            private$flexden$flexden = temp
                            private$flexden$fullden = temp - private$precompute
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
                              temp = pdiscnorm(murep, sd = self$beta, h = self$h) * sum(pi0)
                              ans$d0 = .colSums((temp - private$flexden$flexden) * private$flexden$fullden * self$w, m = self$len, n = length(mu))
                            }
                            if (any(order[2:3] == 1)){
                              temp = dnorm(murep + self$h, sd = self$beta) * private$flexden$fullden * self$w
                            }
                            if (order[2] == 1){
                              ans$d1 = .colSums(temp, m = self$len, n = length(mu)) * -2 * sum(pi0)
                            }
                            if (order[3] == 1){
                              ans$d2 = .colSums(temp * (self$h + murep), m = self$len, n = length(mu)) * -2 * sum(pi0)
                            }

                            ans
                          },
                          computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                            if (any(self$compareattr(mu0, pi0))){
                              self$setflexden(mu0, pi0)
                            }
                            mu0new = c(mu0, newweights)
                            S = cbind(private$flexden$dens, matrix(pdiscnorm(self$data, mean = rep(newweights, rep(self$len, length(newweights))), sd = self$beta, h = self$h), nrow = self$len, ncol = length(newweights)))
                            pi0new = pnnls(S * sqrt(self$w), private$precompute * sqrt(self$w), sum = 1 - sum(self$pi0fixed))$x
                            self$collapsemix(mu0new, pi0new, tol)
                          },
                          estpi0dS = function(){
                            private$S1 = list(a1 = pdiscnorm(self$data, sd = self$beta, h = self$h),
                                              a2 = (2 * cumsum(self$w) - 1) / 2 / sum(self$w))
                          },
                          estpi0d = function(mu0, pi0){
                            ans = vector("list", 2)
                            names(ans) = c("d2", "d3")
                            pn = pnpdiscnorm(self$data, mu0 = mu0, pi0 = pi0, sd = self$beta, h = self$h)
                            S = private$S1$a1 - pn
                            ans$d2 = 2 * sum(S * (pn - private$S1$a2) * self$w); ans$d3 = 2 * sum(S^2 * self$w)
                            ans
                          },
                          estpi0 = function(val = qCvM(0.5, lower.tail = FALSE), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
                            self$modified(pi0 = 1 - tol / 2)
                            self$computemixdist(mix = mix, tol = tol, maxiter = maxiter)
                            r1ll = self$result$ll
                            self$modified(pi0 = 0)
                            self$computemixdist(mix = mix, tol = tol, maxiter = maxiter)
                            nval = sum(self$w) / 3 - sum(self$w * (cumsum(self$w) - 0.5)^2) / sum(self$w)^2

                            if (r1ll + nval < val){
                              self$result = list(iter = 0,
                                       family = self$type,
                                       max.gradient = self$gradientfunction(0, 0, 1, order = c(1, 0, 0))$d0,
                                       mix = list(pt = 0, pr = 1),
                                       ll = self$lossfunction(mu0 = 0, pi0 = 1),
                                       beta = self$beta,
                                       convergence = 0)
                            }else{
                              private$solveestpi0(init = dnpdiscnorm(0, mu0 = self$result$mix$pt, pi0 = self$result$mix$pr, sd = self$beta, h = self$h) / ddiscnorm(0, sd = self$beta, h = self$h),
                                                  val = nval - val, tol = tol, maxiter = maxiter, verbose = verbose)
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
makeobject.npnormcvmw = function(v, mu0, pi0, beta, order = -3){
  npnormcvmw$new(v, mu0, pi0, beta, order = order)
}


