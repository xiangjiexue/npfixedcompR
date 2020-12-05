dnormc = function(x, mean = 0, n, log = FALSE){
  if (any(abs(x) >= 1) | any(abs(mean) >= 1))
    stop("Error in specifying data or mean")
  LLL = max(length(x), length(mean))
  xx = rep(x, length.out = LLL)
  meanx = rep(mean, length.out = LLL)
  dnorm(xx, mean = meanx, sd = (1 - meanx^2) / sqrt(n), log = log)
}

pnormc = function(x, mean = 0, n, lower.tail = TRUE, log.p = FALSE){
  if (any(abs(x) >= 1) | any(abs(mean) >= 1))
    stop("Error in specifying data or mean")
  LLL = max(length(x), length(mean))
  xx = rep(x, length.out = LLL)
  meanx = rep(mean, length.out = LLL)
  pnorm(xx, mean = meanx, sd = (1 - meanx^2) / sqrt(n), lower.tail = lower.tail, log.p = log.p)
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
                        list(pt = seq(from = min(self$data), to = max(self$data), length.out = 10),
                             pr = rep(1 / 10, 10))
                      },
                      beta = NULL,
                      type = "npnormc",
                      initialize = function(data, mu0, pi0, beta){
                        self$data = data
                        self$len = length(data)
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$flexden = list(flexden = numeric(self$len), fullden = numeric(self$len))
                        attr(private$flexden, "mu0") = NULL
                        attr(private$flexden, "pi0") = NULL
                        private$precompute = dnpnormc(data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, n = self$beta)
                      },
                      modified = function(mu0, pi0, beta){
                        if (!missing(mu0)) self$mu0fixed = mu0
                        if (!missing(pi0)) self$pi0fixed = pi0
                        if (!missing(beta)) self$beta = beta
                        private$precompute = dnpnormc(data, mu0 = self$mu0fixed, pi0 = self$pi0fixed, n = self$beta)
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
                          private$flexden$dens = dnormc(self$data, mean = rep(mu0, rep(self$len, length(mu0))), n = self$beta)
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
                        ans = vector("list", 3)
                        names(ans) = c("d0", "d1", "d2")
                        if (length(mu) > 0){
                          murep = rep(mu, rep(self$len, length(mu)))
                          temp = dnormc(self$data, mean = murep, n = self$beta)
                          private$flexden$temp = matrix(temp, ncol = length(mu))
                          # only d0 is implemented
                          if (order[1] == 1){
                            ans$d0 = .colSums((private$flexden$flexden - temp * sum(pi0)) / private$flexden$fullden, m = self$len, n = length(mu))
                          }
                          if (order[2] == 1){
                            temp2 = ((self$beta + 4) * murep^3 - 3 * self$beta * murep^2 * self$data + murep * (2 * self$beta * self$data^2 + self$beta - 2) -
                                       self$beta * self$data - 2 * murep^5) / (1 - murep^2)^3
                            ans$d1 = .colSums(temp2 * temp / private$flexden$fullden, m = self$len, n = length(mu))  * sum(pi0)
                          }

                        }

                        ans
                      },
                      computeweights = function(mu0, pi0, newweights, tol = 1e-6){
                        if (any(self$compareattr(mu0, pi0))){
                          self$setflexden(mu0, pi0)
                        }
                        if (length(newweights) > 0){
                          mu0new = c(mu0, newweights)
                          pi0new = c(pi0, rep(0, length(newweights)))
                          sp = cbind(private$flexden$dens, private$flexden$temp)
                        }else{
                          mu0new = mu0
                          pi0new = pi0
                          sp = private$flexden$dens
                        }
                        fp = private$flexden$fullden
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
                      mflag = "d1",
                      S1 = NULL
                    ))

#' @rdname makeobject
#' @export
makeobject.npnormcll = function(v, mu0, pi0, beta, order){
  npnormcll$new(v, mu0, pi0, beta)
}
