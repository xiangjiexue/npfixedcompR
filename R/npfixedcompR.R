#' npfixedcompR
#'
#' npfixedcompR
#'
#' @docType package
#' @author Xiangjie Xue
#' @importFrom stats dnorm uniroot pnorm dt optimise pt cov cov2cor runif
#' @importFrom grDevices rainbow
#' @importFrom graphics hist rect lines plot points abline
#' @importFrom lsei pnnls pnnqp
#' @importFrom goftest qCvM qAD
#' @importFrom nloptr nloptr
#' @importFrom R6 R6Class
#' @name npfixedcompR
NULL

#' These functions are used to make the object for computing the non-paramtric mixing
#' distribution or estimating the proportion of zero using non-parametric methods.
#'
#' This is a S3 generic function for making the object for computing the non-parametric
#' mixing distribution or estimating the proportion of zero.
#'
#' currently implemented families are:
#'
#' - npnormll : normal density using maximum likelihood
#'
#' - npnormllw : Binned version of normal density using maximum likelihood
#'
#' - npnormcvm : normal density using Cramer-von Mises
#'
#' - npnormcvmw : Binned version of normal density using Cramer-von Mises
#'
#' - npnormad : normal density using Anderson-Darling
#'
#' - npnormadw : Binned version of normal density using Anderson-Darling
#'
#' - nptll : t density using maximum likelihood
#'
#' - npnormcll : one-parameter normal distribution used for approximate sample
#' correlation coefficients using maximum likelihood. This does not have a
#' corresponding estimation of zero due to incomplete of theory.
#'
#' The default method used is npnormll.
#'
#' @title Making object for computation
#' @param v the object either numeric or the implmented family
#' @param mu0 A vector of support points
#' @param pi0 A vector of weights corresponding to the support points
#' @param beta structual parameter.
#' @param order the parameter for the binned version.
#' @param method An implemented family; see details
#' @param ... other parameter passed to the constructor.
#' @rdname makeobject
#' @export
makeobject = function(v, method = "npnormll", ...){
  f = match.fun(paste0("makeobject.", method))
  f(v = v, ...)
}

# Taken from nspmix with modification
# For simplifying the mixing distribution
unique.disc = function(pt, pr, prec=0) {
  k = sum(pr)
  if (length(pt) == 1 ) return(list(pt = pt, pr = pr))
  prec = rep(prec, len=2)
  j0  = pr <= prec[2]
  pt = pt[!j0]
  pr = pr[!j0]
  repeat {
    i = diff(pt) <= prec[1]                # will use all odd indexes of TRUEs
    if(sum(i) == 0) break
    noti = c(!i, TRUE)                              # not i
    even = seq(2, by=2, len=length(i))
    i[even] = i[even] & noti[even-1] & noti[even+1] # use evens if only isolated
    ind = which(i)
    pt[ind] = (pt[ind] * pr[ind] + pt[ind+1] * pr[ind+1]) /
      (pr2 <- pr[ind] + pr[ind+1])                  # collapse neighbouring pairs
    pr[ind] = pr2
    pt = pt[-(ind+1)]                               # remove the second of a pair
    pr = pr[-(ind+1)]
  }
  index = order(pt, decreasing = FALSE)
  list(pt = pt[index], pr = pr[index] / sum(pr) * k)
}

# The following two functions taken from nspmix with modification
whist = function(x, w=1, breaks="Sturges", plot=TRUE, freq=NULL,
                 xlim=NULL, ylim=NULL, xlab="Data", ylab=NULL, main=NULL,
                 add=FALSE, col=NULL, border=NULL, lwd=1, ...) {
  r = hist(x, breaks=breaks, plot=FALSE)
  breaks = r$breaks
  i = cut(x, breaks, include.lowest=TRUE)
  f = tapply(rep(w,len=length(i)), i, sum)            # frequency
  f[is.na(f)] = 0
  dimnames(f)[[1]] = NULL
  d = f / sum(f) / (breaks[2] - breaks[1])            # density
  if(! is.null(freq) && ! freq) {
    y = d
    if(is.null(ylab)) ylab = "Density"
  }
  else {
    y = f
    if(is.null(ylab)) ylab = "Frequency"
  }
  ny = length(y)
  if (is.null(xlim)) xlim = range(breaks)
  if(is.null(ylim)) ylim = range(0, y, finite=TRUE)
  else {
    ymax = max(y)
    ylim = c(0, pmin(pmax(ymax, max(ylim)), 2*ymax))
  }
  if(plot) {
    if(!add) plot(r$mids, y, xlim=xlim, ylim=ylim, type="n", frame.plot=FALSE,
                  xlab=xlab, ylab=ylab, main=main, ...)
    rect(breaks[-(ny+1)], 0, breaks[-1], y, col=col, border=border, lwd=lwd)
    lines(range(breaks), c(0,0), col=border)
  }
  else list(breaks=breaks, counts=f, density=d,
            mids=breaks[-(ny+1)] + diff(breaks) * .5)
}

#' Bin the continuous data.
#'
#' This function bins the continuous data using the equal-width bin in order to
#' speed up some functions in this package.
#'
#' h is taken as 10^(-order)
#'
#' the observations are rounded down to the bins ..., -h, 0, h, ...
#'
#' To further speed up the process, omit the bin that has 0 count.
#'
#' @title Bin the continuous data.
#' @param data the observation to be binned.
#' @param order see the details
#' @return a list with v be the representative value of each bin and w be the count in the corresponding bin.
#' @export
bin = function(data, order = -5){
  h = 10^order
  data = floor(data / h) * h
  t = table(data)
  index = t != 0
  list(v = as.numeric(names(t))[index], w = as.numeric(t)[index])
}

npfixedcompR = R6::R6Class("npfixedcompR",
                       public = list(
                         mu0fixed = 0,
                         pi0fixed = 0,
                         data = NULL,
                         len = NULL,
                         printmethodflag = function(){
                           private$methodflag
                         },
                         printgridpoints = function(){
                           if (is.null(private$gridpoints)){
                             if (is.null(x$setgridpoints)){
                               private$gridpoints = seq(from = min(self$data), to = max(self$data), length = 100)
                             }else{
                               self$setgridpoints()
                             }
                           }
                           private$gridpoints
                         },
                         checklossfunction = function(mu0, pi0, eta, p, maxit = 100){
                           llorigin = self$lossfunction(mu0, pi0)
                           con = sum(p * eta)
                           sigma = 0.5; alpha = 1/3; u = -1;
                           repeat{
                             u = u + 1
                             pi0new = pi0 + sigma^u * eta
                             lhs = -self$lossfunction(mu0, pi0new)
                             rhs = -llorigin + alpha * sigma^u * con
                             if (lhs > rhs){
                               pi0 = pi0new
                               break
                             }
                             if (u + 1 > maxit)
                               break
                           }
                           index = order(mu0, decreasing = FALSE)
                           unique.disc(mu0[index], pi0[index])
                         },
                         collapsemix = function(mu0, pi0, tol){
                           ll = self$lossfunction(mu0, pi0)
                           ntol = max(tol * 0.1, ll * 1e-16)
                           repeat{
                             if (length(mu0) <= 1)
                               break
                             prec = c(10 * min(diff(mu0)))
                             r = unique.disc(mu0, pi0, c(prec, 0))
                             nll = self$lossfunction(r$pt, r$pr)
                             if (nll <= ll + ntol){
                               mu0 = r$pt; pi0 = r$pr
                             }else{
                               break
                             }
                           }
                           list(pt = mu0, pr = pi0)
                         },
                         compareattr = function(mu0, pi0){
                           # move this to general constructor
                           if (is.null(private$flexden)){
                             TRUE
                           }else{
                             if (length(mu0) != length(attr(private$flexden, "mu0")) | length(pi0) != length(attr(private$flexden, "pi0"))){
                               TRUE
                             }else{
                               sum((mu0 - attr(private$flexden, "mu0"))^2) > 1e-12 | sum((pi0 - attr(private$flexden, "pi0"))^2) > 1e-12
                             }
                           }
                         },
                         estpi0d = function(mu0, pi0){
                           ans = vector("list", 2)
                           names(ans) = c("d2", "d3")
                           ans$d2 = x$gradientfunction(mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                           ans
                         }
                       ),
                       private = list(
                         gridpoints = NULL
                       ))

#' This is a S3 generic function for computing non-parametric mixing distribution
#'
#' The full list of implemented family is in \code{\link{makeobject}}.
#'
#' @title Computing non-parametric mixing distribution
#' @param x a object from implemented family generated by \code{\link{makeobject}}.
#' @param mix initial mixing distribution (proper)
#' @param tol tolerance to stop the code
#' @param maxiter maximum iteration to stop the code
#' @param verbose logical. If TRUE, the intermediate results will be printed.
#' @examples
#' data = rnorm(500, c(0, 2))
#' x = makeobject(data, method = "npnormll")
#' computemixdist(x)
#' x = makeobject(data, method = "npnormcvm")
#' computemixdist(x)
#' x = makeobject(data, method = "npnormad")
#' computemixdist(x)
#' @export
computemixdist = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  UseMethod("computemixdist")
}

#' @export
computemixdist.npfixedcompR = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    mix = x$initpoints()
  }
  mu0 = mix$pt; pi0 = mix$pr
  pi0 = pi0 * (1 - sum(x$pi0fixed))
  iter = 0; convergence = 0
  closs = x$lossfunction(mu0, pi0)

  repeat{
    r = x$computeweights(mu0, pi0, solvegradientmultiple(x, mu0, pi0), tol = tol)
    mu0 = r$pt; pi0 = r$pr
    iter = iter + 1
    nloss = x$lossfunction(mu0, pi0)

    if (verbose){
      cat("Iteration: ", iter, "\n")
      cat(paste0("Support Point ", round(r$pt, -ceiling(log10(tol))), " with probability ", round(r$pr, -ceiling(log10(tol))), "\n"))
      cat("Current log-likelihood ", as.character(round(-nloss, -ceiling(log10(tol)))), "\n")
    }

    if (closs - nloss < tol){
      convergence = 0
      break
    }

    if (iter > maxiter){
      convergence = 1
      break
    }
    closs = nloss
  }

  mu0new = c(mu0, x$mu0fixed); pi0new = c(pi0, x$pi0fixed)
  index = order(mu0new, decreasing = FALSE)
  r = unique.disc(mu0new[index], pi0new[index])

  ans = list(iter = iter,
             family = x$type,
             max.gradient = min(x$gradientfunction(mu0, mu0, pi0, order = c(1, 0, 0))$d0),
             mix = list(pt = r$pt, pr = r$pr),
             ll = nloss,
             beta = x$beta,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

#' computing non-parametric mixing distribution with estimated proportion at 0
#'
#' This is a S3 generic function for computing non-parametric mixing
#' distribution with estimated proportion at 0. Different families will
#' have different threshold values.
#'
#' The full list of implemented family is in \code{\link{makeobject}}.
#'
#' @title Computing non-parametric mixing distribution with estimated proportion at 0
#' @param x a object from implemented family
#' @param val the threshold value.
#' @param mix initial mixing distribution (proper)
#' @param tol tolerance to stop the code
#' @param maxiter maximum iteration to stop the code
#' @param verbose logical. If TRUE, the intermediate results will be printed.
#' @param ... parameters passed to the specific method.
#' @examples
#' data = rnorm(500, c(0, 2))
#' pi0 = 0.5
#' x = makeobject(data, method = "npnormll")
#' estpi0(x)
#' x = makeobject(data, method = "npnormcvm")
#' estpi0(x)
#' x = makeobject(data, method = "npnormad")
#' estpi0(x)
#' @export
estpi0 = function(x, ...){
  f = match.fun(paste0("estpi0.", class(x)[1]))
  f(x = x, ...)
}
