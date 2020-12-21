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

#' This documentation gives a detailed description of the npfixedcompR object. This implementation uses
#' the R6 object.
#'
#' The structure of the object currently has three layers: The foundation layer is the
#' npfixedcompR class, which has components common to all the specific classes; The second
#' layer is the distributional layer, which has components common to all the classes with
#' the same distribution (the npnorm class); The third layer is the specific layer
#' (the npnormll class), which contains components specific to this particular loss. Since
#' the nptll and the npnormcll class only has one loss implemented, they are implemented in
#' the second layer.
#'
#' For the following, if a component is listed as public, then it can be accessed as needed.
#' If a component is listed as private, then it can not be accessed. There might be functions
#' listed as public can be used to modified the private values.
#'
#' The component in the npfixedcompR class are as follows.
#'
#' - \code{mu0fixed} (public) : The vector of the locations of fixed support points.
#' This can be set by \code{initialize} or \code{modified} implemented in the last year.
#'
#' - \code{pi0fixed} (public) : The vector of the weights of fixed support points.
#' This can be set by \code{initialize} or \code{modified} implemented in the last year.
#'
#' - \code{data} (public) : The observations. The observations can only be set by
#' \code{\link{makeobject}} via \code{initialzed}. Once the object is defined, the observations
#' can not be modified.
#'
#' - \code{len} (public) : The length of the observations. This is different to the effective
#' sample size when dealing with the binned version.
#'
#' - \code{result} (public) : The estimated mixing distribution. This is used to store the result
#' computed by \code{computemixdist} or \code{estpi0}.
#'
#' - \code{methodflag} (public) : This function is used to print or set the algorithm used for
#' finding the new support points. The default argument is \code{NULL}, which prints the private
#' component \code{mflag} implemented in the final layer. Other inputs can be \code{d0}, \code{d1}
#' and \code{d2}, which change the algorithm used. The algorithm used should depend on the final layer
#' (hence \code{mflag} is implemented in the final layer).
#'
#' - \code{checklossfunction} (public) : The function used for Armijo rule; see Wang (2007).
#'
#' - \code{collapsemix} (public) : The function used for collapsing mixing distributions; see \code{\link[nspmix]{cnm}}.
#'
#' - \code{compareattr} (public) : The function used for comparing the attributes \code{mu0} and \code{pi0} for performence
#' purposes. If the geometry of \code{mu0} and \code{pi0} is close to the stored ones, there is no need to recompute
#' the information related to this pair.
#'
#' - \code{solvegradientmultiple} (public) : The function used for computing the new support points.
#' This is function calls for \code{solvegradientmultipled0}, \code{solvegradientmultipled1} and
#' \code{solvegradientmultipled2} according to the private component \code{mflag}. The detailed descriptions
#' of the algorithms used can be found in Appendix B of the thesis. The break points of the smaller intervals
#' are computed through \code{printgridpoints} implemented in the distributional layer.
#'
#' - \code{computemixdist} (public) : The function used for computing the mixing distribution given
#' \code{mu0fixed} and \code{pi0fixd}. The result is stored in \code{result} component.
#'
#' - \code{estpi0d} (public) : The function used for computing the derivative when estimate the null proportion.
#' This default implementation can be used in the Newton-Raphson method. This is overwritten by the
#' espti0d in the last layer, which can be used in the Halley method.
#'
#' - \code{solvegradientmultipled0} (private), \code{solvegradientmultipled1} (private) and
#' \code{solvegradientmultipled2} (private) : The functions for computing the new support points. The underlying
#' computing functions \code{solvegradientsingled1} and \code{solvegradientsingled2} are vectorised
#' implementations while \code{solvegradientsingled0} uses the \code{\link[stats]{optimize}}.
#'
#' - \code{solveestpi0} (private) : The function for computing the null proportion. The algofithm used depends
#' on \code{estpi0d}. This is essentially the root-finding algorithm and the formulation can be found in
#' Appendix B of the thesis.
#'
#' - \code{gridpoints} (private) : The break points for the smaller interval when computing the new support
#' points. This is initialised by \code{printgridpoints} in the distributional layer.
#'
#' The components in the distributional layer are as follows.
#'
#' - \code{setgridpoints} (public) : The function used to initialise \code{gridpoints}. The range of the
#' support can be found in Chapter 3 of the thesis.
#'
#' - \code{initpoints} (public) : The function for generating the starting mixing distribution, if not
#' specified in \code{computemixdist}; see \code{link[nspmix]{npnorm}} for example.
#'
#' The componets in the final layer are as follows.
#'
#' - \code{beta} (public) : The structural parameter. The structural parameter is always considered as fixed
#' and can be changed using \code{initisalize} and \code{modified}.
#'
#' - \code{type} (public) : The flag for component distributions.
#'
#' - \code{initialize} (public) : The function for initialising the object for computing. It sets the data,
#' the length, fixed support points, structural parameter and any precomputed values.
#'
#' - \code{modified} (public) : The function for modifying the object for further computing. It sets the
#' fixed support points, the structural parameter and any precomputed values related to the change of
#' the fixed support points.
#'
#' - \code{lossfunction} (public) : THe function for computing the loss function given a mixing distribution.
#'
#' - \code{setflexden} (public) : This function is used to set the precomputed values for performance
#' purposes.
#'
#' - \code{gradientfunction} (public) : This function computes the gradient for the implemented families.
#' The output is a list of length 3: The gradient \code{d0}, the first derivative with respect to the location
#' parameter \code{d1} and the second derivateive with respect to the location parameter \code{d2}.
#'
#' - \code{computeweights} (public) : This function computes the new weights given fixed support points,
#' which is needed in \code{computemixdist} in the first layer.
#'
#' - \code{estpi0dS} (public) : Generating precomputed values for estimating the null proportion. This is
#' primarily for performance purposes.
#'
#' - \code{estpi0d} (public) : This function overwrites the function of the same name in the first layer
#' for a faster computation of the null proportions.
#'
#' - \code{estpi0} (public) : This function computes the null proportions depending on the given loss/distance.
#'
#' - \code{precompute} (private), \code{flexden} (private), \code{S1} (private) : The object for storing the
#' precomputed values for performance purposes.
#'
#' @title The npfixedcompR object
#' @references
#' Wang, Yong. "On Fast Computation of the Non-Parametric Maximum Likelihood Estimate of a Mixing Distribution." Journal of the Royal Statistical Society. Series B (Statistical Methodology) 69, no. 2 (2007): 185-98. \url{http://www.jstor.org/stable/4623262}.
#' @name npfixedcompRobject
NULL

#' These functions are used to make the object for computing the non-paramtric mixing
#' distribution or estimating the proportion of zero using non-parametric methods.
#'
#' This is a generic function for making the object for computing the non-parametric
#' mixing distribution or estimating the proportion of zero.
#'
#' current implemented families are:
#'
#' - npnormll : normal density using maximum likelihood (Chapter 3). The default beta is 1.
#'
#' - npnormllw : Binned version of normal density using maximum likelihood (Chapter 6).
#' The default beta is 1 and the default order is -3.
#'
#' - npnormcvm : normal density using the Cramer-von Mises distance (Chapter 5). The default beta is 1.
#'
#' - npnormcvmw : Binned version of normal density using the Cramer-von Mises distance (Chapter 6).
#' The default beta is 1 and the default order is -3.
#'
#' - npnormad : normal density using the Anderson-Darling distance (Chapter 5). The default beta is 1.
#'
#' - npnormadw : Binned version of normal density using the Anderson-Darling distance (Chapter 6)
#' The default beta is 1 and the default order is -3.
#'
#' - nptll : t-density using maximum likelihood (Chapter 3). The default beta is infinity (normal distribution).
#'
#' - npnormcll : the one-parameter normal distribution used for approximating the sample
#' correlation coefficients using maximum likelihood. This does not have a
#' corresponding estimation of zero due to incompleted theory (Chapter 8).
#' There is no default beta. The structure beta is the number of observations.
#'
#' The default method used is npnormll.
#'
#' The detailed description of the npfixedcompRobject class is in \code{\link{npfixedcompRobject}}.
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
  breaks = pmin(pmax(r$breaks, min(x)), max(x)) # the breaks can not exceed the range
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
                         result = NULL,
                         methodflag = function(methodflag = NULL){
                           if (!is.null(methodflag)) private$mflag = methodflag
                           private$mflag
                         },
                         printgridpoints = function(){
                           if (is.null(private$gridpoints)){
                             if (is.null(self$setgridpoints)){
                               private$gridpoints = seq(from = min(self$data), to = max(self$data), length = 100)
                             }else{
                               self$setgridpoints()
                             }
                           }
                           private$gridpoints
                         },
                         checklossfunction = function(mu0, pi0, eta, p, maxit = 100){
                           llorigin = self$lossfunction(mu0, pi0)
                           sigma = 0.5; alpha = 1/3; u = -1;
                           con = -sum(p * eta)
                           repeat{
                             u = u + 1
                             pi0new = pi0 + sigma^u * eta
                             lhs = self$lossfunction(mu0, pi0new)
                             rhs = llorigin + alpha * sigma^u * con
                             if (lhs < rhs){
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
                         solvegradientmultiple = function(mu0, pi0, tol = 1e-6){
                           points = self$printgridpoints()
                           switch(private$mflag,
                                  "d1" = private$solvegradientmultipled1(mu0, pi0, tol),
                                  "d2" = private$solvegradientmultipled2(mu0, pi0, tol),
                                  "d0" = private$solvegradientmultipled0(mu0, pi0, tol))
                         },
                         computemixdist = function(mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
                           if (is.null(mix)){
                             mix = self$initpoints()
                           }
                           mu0 = mix$pt; pi0 = mix$pr
                           pi0 = pi0 * (1 - sum(self$pi0fixed))
                           iter = 0; convergence = 0
                           closs = self$lossfunction(mu0, pi0)

                           repeat{
                             newpoints = self$solvegradientmultiple(mu0, pi0)
                             r = self$computeweights(mu0, pi0, newpoints, tol = tol)
                             mu0 = r$pt; pi0 = r$pr
                             iter = iter + 1
                             nloss = self$lossfunction(mu0, pi0)

                             if (verbose){
                               cat("Iteration: ", iter, "\n")
                               cat(paste0("Support Point ", round(r$pt, -ceiling(log10(tol))), " with probability ", round(r$pr, -ceiling(log10(tol))), "\n"))
                               cat("Current Loss ", as.character(round(nloss, -ceiling(log10(tol)))), "\n")
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

                           mu0new = c(mu0, self$mu0fixed); pi0new = c(pi0, self$pi0fixed)
                           index = order(mu0new, decreasing = FALSE)
                           r = unique.disc(mu0new[index], pi0new[index])

                           self$result = list(iter = iter,
                                      family = self$type,
                                      max.gradient = min(self$gradientfunction(mu0, mu0, pi0, order = c(1, 0, 0))$d0),
                                      mix = list(pt = r$pt, pr = r$pr),
                                      ll = nloss,
                                      beta = self$beta,
                                      convergence = convergence)

                           attr(self$result, "class") = "nspmix"
                         },
                         estpi0d = function(mu0, pi0){
                           ans = vector("list", 2)
                           names(ans) = c("d2", "d3")
                           ans$d2 = x$gradientfunction(mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                           ans
                         }
                       ),
                       private = list(
                         gridpoints = NULL,
                         solvegradientsingled0 = function(mu0, pi0, lower, upper, tol = 1e-6){
                           optimise(function(mu_){
                             self$gradientfunction(mu = mu_, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                           }, lower = lower, upper = upper, tol = tol)$minimum
                         },
                         solvegradientmultipled0 = function(mu0, pi0, tol = 1e-6){
                           points = private$gridpoints
                           pointsval = diff(self$gradientfunction(mu = points, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0)
                           index = seq(1, by = 1, length = length(points) - 2)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
                           if (length(index) >= 1){
                             r = sapply(index, function(ddd){
                               private$solvegradientsingled0(mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 2], tol)
                             })
                             r = c(points[1], r, points[length(points)])
                             dd = self$gradientfunction(mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                             r = r[dd < 0]
                           }else{
                             dd = self$gradientfunction(mu = c(points[1], points[length(points)]), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                             r = c(points[1], points[length(points)])[dd < 0]
                           }
                           private$flexden$temp = private$flexden$temp[, dd < 0, drop = FALSE]
                           r
                         },
                         solvegradientsingled1 = function(mu0, pi0, lower, upper, tol = 1e-6){
                           # vectorised implementation
                           a = lower; b = upper
                           f1 = self$gradientfunction(mu = c(a, b), mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
                           fa = f1[1:length(a)]; fb = f1[(length(a) + 1) : length(f1)]

                           fs = fb; s = b
                           fc = numeric(length(a))
                           index = rep(FALSE, length(a))
                           repeat{
                             if (sum(index) == length(a))
                               break

                             index = abs(fb) < tol | abs(fs) < tol | abs(b - a) < tol

                             c1 = (a + b) / 2
                             fc[!index] = self$gradientfunction(c1[!index], mu0, pi0, order = c(0, 1, 0))$d1

                             # does not matter how the following changes for a[index], b[index], etc. changes, not referenced as long as
                             # s is not assigned to a or b
                             j0 = fa != fc & fb != fc
                             s[j0] = a[j0] * fb[j0] * fc[j0] / (fa[j0] - fb[j0]) / (fa[j0] - fc[j0]) +
                               b[j0] * fa[j0] * fc[j0] / (fb[j0] - fa[j0]) / (fb[j0] - fc[j0]) +
                               c1[j0] * fa[j0] * fb[j0] / (fc[j0] - fa[j0]) / (fc[j0] - fb[j0])
                             s[!j0] = b[!j0] - fb[!j0] * (b[!j0] - a[!j0]) / (fb[!j0] - fa[!j0])

                             fs[!index] = self$gradientfunction(s[!index], mu0, pi0, order = c(0, 1, 0))$d1

                             j0 = fs > 0
                             j1 = s > a
                             a[!j0 & !index & j1] = s[!j0 & !index & j1]
                             fa[!j0 & !index & j1] = fs[!j0 & !index & j1]
                             j1 = s < b
                             b[j0 & !index & j1] = s[j0 & !index & j1]
                             fb[j0 & !index & j1] = fs[j0 & !index & j1]

                             j0 = fc > 0
                             a[!j0 & !index] = c1[!j0 & !index]
                             fa[!j0 & !index] = fc[!j0 & !index]
                             b[j0 & !index] = c1[j0 & !index]
                             fb[j0 & !index] = fc[j0 & !index]
                           }

                           b
                         },
                         solvegradientmultipled1 = function(mu0, pi0, tol = 1e-6){
                           points = private$gridpoints
                           pointsval = self$gradientfunction(mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
                           index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
                           if (length(index) >= 1){
                             # r = sapply(index, function(ddd){
                             #   solvegradientsingled1(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
                             # })
                             r = c(points[1], private$solvegradientsingled1(mu0 = mu0, pi0 = pi0, lower = points[index], upper = points[index + 1], tol), points[length(points)])
                             dd = self$gradientfunction(mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                             r = r[dd < 0]
                           }else{
                             dd = self$gradientfunction(mu = c(points[1], points[length(points)]), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                             r = c(points[1], points[length(points)])[dd < 0]
                           }
                           private$flexden$temp = private$flexden$temp[, dd < 0, drop = FALSE]
                           r
                         },
                         solvegradientsingled2 = function(mu0, pi0, lower, upper, tol = 1e-6){
                           # This implements the vectorised version
                           neg = lower; pos = upper; init = (neg + pos) / 2
                           r = self$gradientfunction(init, mu0, pi0, order = c(0, 1, 1))
                           d1 = r$d1
                           d2 = r$d2
                           d = numeric(length(lower))
                           mflag = rep(FALSE, length(neg)) # There is no need to calculate the gradient function for mid-point again.

                           index = rep(FALSE, length(neg))
                           repeat{
                             if (sum(index) == length(neg))
                               break
                             index = abs(d1) < tol | pos - neg < tol

                             if (sum(mflag) > 0){
                               s = (pos + neg) / 2
                               d[!index & mflag] = self$gradientfunction(s[!index & mflag], mu0, pi0, order = c(0, 1, 0))$d1

                               j0 = d < 0
                               neg[j0 & !index & mflag] = s[j0 & !index & mflag]
                               pos[!j0 & !index & mflag] = s[!j0 & !index & mflag]
                             }

                             j0 = d1 < 0
                             neg[j0] = pmax(neg[j0], init[j0])
                             pos[!j0] = pmin(pos[!j0], init[!j0])

                             init[!index] = init[!index] - d1[!index] / d2[!index]
                             mflag = init > neg & init < pos
                             init[!index & !mflag] = (pos[!index & !mflag] + neg[!index & !mflag]) / 2

                             r = self$gradientfunction(init[!index], mu0, pi0, order = c(0, 1, 1))
                             d1[!index] = r$d1
                             d2[!index] = r$d2
                           }

                           init
                         },
                         solvegradientmultipled2 = function(mu0, pi0, tol = 1e-6){
                           points = private$gridpoints
                           pointsval = self$gradientfunction(mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
                           index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
                           if (length(index) >= 1){
                             r = c(points[1], private$solvegradientsingled2(mu0 = mu0, pi0 = pi0, lower = points[index], upper = points[index + 1], tol), points[length(points)])
                             dd = self$gradientfunction(mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                             r = r[dd < 0]
                           }else{
                             dd = self$gradientfunction(mu = c(points[1], points[length(points)]), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
                             r = c(points[1], points[length(points)])[dd < 0]
                           }
                           private$flexden$temp = private$flexden$temp[, dd < 0, drop = FALSE]
                           r
                         },
                         solveestpi0 = function(init, val, tol = 1e-6, maxiter = 100, verbose = FALSE){
                           neg = 0; pos = 1 - tol / 2
                           iter = 1;
                           self$modified(pi0 = init)
                           self$computemixdist(mix = self$result$mix, maxiter = maxiter, tol = tol)
                           d1 = self$result$ll + val
                           d = self$estpi0d(self$result$mix$pt, self$result$mix$pr)
                           d2 = d$d2; d3 = d$d3
                           repeat{
                             if (verbose){
                               cat("Iteration", iter, "\n")
                               cat("Current pi0 value =", round(init, -ceiling(log10(tol))), "with function value", d1,
                                   "(iterations", self$result$iter, ")\n")
                               cat("Current Bracket: lower =", neg, "upper =", pos, "\n")
                             }

                             if (iter > maxiter | abs(d1) < tol | pos - neg < tol)
                               break
                             iter = iter + 1

                             if (d1 < 0){
                               if ((init - neg) / (pos - neg) < 0.5){
                                 s = (pos + neg) / 2
                                 self$modified(pi0 = s)
                                 self$computemixdist(mix = self$result$mix, maxiter = maxiter, tol = tol)

                                 if (self$result$ll + val < 0){
                                   neg = s
                                 }else{
                                   pos = s
                                 }
                               }
                               neg = max(neg, init)
                             }else{
                               if ((pos - init) / (pos - neg) < 0.5){
                                 s = (pos + neg) / 2
                                 self$modified(pi0 = s)
                                 self$computemixdist(mix = self$result$mix, maxiter = maxiter, tol = tol)

                                 if (self$result$ll + val < 0){
                                   neg = s
                                 }else{
                                   pos = s
                                 }
                               }
                               pos = min(pos, init)
                             }

                             if (is.null(d3)){
                               init = init - d1 * (1 - init) / d2;
                             }else{
                               init = init - 2 * d1 * d2 * (1 - init) / (2 * d2^2 - d1 * d3)
                             }


                             if (init < neg | init > pos){
                               init = (pos + neg) / 2
                             }

                             self$modified(pi0 = init)
                             self$computemixdist(mix = self$result$mix, maxiter = maxiter, tol = tol)
                             d1 = self$result$ll + val
                             d = self$estpi0d(self$result$mix$pt, self$result$mix$pr)
                             d2 = d$d2; d3 = d$d3
                           }
                         }
                       ))

#' computing non-parametric mixing distribution with estimated proportion at 0
#'
#' This is a function for computing non-parametric mixing
#' distribution with estimated proportion at 0. Different families will
#' have different threshold values.
#'
#' The parameters are listed as follows:
#'
#' - val: Threshold value
#'
#' - mix: The initial proper mixing distribution.
#'
#' - tol: tolerance to stop the code.
#'
#' - maxiter: maximum iterations allowed.
#'
#' - verbose: logical; Whether to print the intermediate results.
#'
#' It is not shown in the parameter section since various method have different
#' default thresold values and this function essentially calls the class method
#' in the object.
#'
#' The full list of implemented families is in \code{\link{makeobject}}.
#'
#' @title Computing non-parametric mixing distribution with estimated proportion at 0
#' @param x a object from implemented family
#' @param ... parameters above passed to the specific method.
#' @examples
#' data = rnorm(500, c(0, 2))
#' pi0 = 0.5
#' x = makeobject(data, method = "npnormll")
#' estpi0(x)
#' x = makeobject(data, method = "npnormllw")
#' estpi0(x)
#' x = makeobject(data, method = "npnormcvm")
#' estpi0(x)
#' x = makeobject(data, method = "npnormcvm")
#' estpi0(x)
#' x = makeobject(data, method = "npnormad")
#' estpi0(x)
#' x = makeobject(data, method = "nptll")
#' estpi0(x)
#' @export
estpi0 = function(x, ...){
  x$estpi0(...)
  x$result
}

#' Computing non-parametric mixing distribution
#'
#' The full list of implemented family is in \code{\link{makeobject}}.
#'
#' The avaliable parameters are listed as follows:
#'
#' - mix: The initial proper mixing distribution
#'
#' - tol: tolerance to stop the code
#'
#' - maxiter: maximum iterations allowed.
#'
#' - verbose: logical; whether to print the intermediate results.
#'
#' This function essentially calls the class method in the object.
#'
#' @title Computing non-parametric mixing distribution
#' @param x a object from implemented family generated by \code{\link{makeobject}}.
#' @param ... parameters above passed to the specific method
#' @examples
#' data = rnorm(500, c(0, 2))
#' pi0 = 0.5
#' x = makeobject(data, pi0 = pi0, method = "npnormll")
#' computemixdist(x)
#' x = makeobject(data, pi0 = pi0, method = "npnormllw")
#' computemixdist(x)
#' x = makeobject(data, pi0 = pi0, method = "npnormcvm")
#' computemixdist(x)
#' x = makeobject(data, pi0 = pi0, method = "npnormcvmw")
#' computemixdist(x)
#' x = makeobject(data, pi0 = pi0, method = "npnormad")
#' computemixdist(x)
#' x = makeobject(data, pi0 = pi0, method = "nptll")
#' computemixdist(x)
#' @export
computemixdist = function(x, ...){
  x$computemixdist(...)
  x$result
}
