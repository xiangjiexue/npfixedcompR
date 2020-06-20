#' npfixedcompR
#'
#' npfixedcompR
#'
#' @docType package
#' @author Xiangjie Xue
#' @importFrom stats dnorm uniroot pnorm
#' @importFrom graphics hist
#' @importFrom lsei pnnls
#' @importFrom goftest qCvM
#' @name npfixedcompR
NULL

#' make object for computation
#'
#' This is a S3 generic function for making the object for computing the non-parametric
#' mixing distribution
#'
#' The default method used is npnormfcll.
#'
#' @title make object for computation
#' @param v the object either numeric or the implmented family
#' @param mu0 A vector of support points
#' @param pi0 A vector of weights corresponding to the support points
#' @param beta structual parameter
#' @export
makeobject = function(v, mu0, pi0, beta){
  UseMethod("makeobject", v)
}

#' @rdname makeobject
makeobject.default = function(v, mu0 = 0, pi0 = 0, beta = 1){
  makeobject.npnormll(v, mu0, pi0, beta)
}

#' loss function
#'
#' This is a S3 generic loss function for computing minimisation problem in order to find
#' non-parametric estimator
#'
#' @title loss function
#' @param x a object from implemented family
#' @param mu0 a vector of support points
#' @param pi0 a vector of weights corresponding to the support points
#' @export
lossfunction = function(x, mu0, pi0){
  UseMethod("lossfunction", x)
}


#' gradient function
#'
#' This is a S3 generic function returning gradident function with order specified as a list
#'
#' @title The gradient function
#' @param x a object from implemented family
#' @param mu a vector of locations to look for
#' @param mu0 a vector of support points
#' @param pi0 a vector of weights corresponding to the support points
#' @param order a vector of length 3 with either 1 or 0. 1 to compute, 0 not to compute.
#' @return a list of length 3. d0, d1, d2 correspond to the order.
#' @export
gradientfunction = function(x, mu, mu0, pi0, order){
  UseMethod("gradientfunction", x)
}

#' computing non-parametric mixing distribution
#'
#' This is a S3 generic function for computing non-parametric mixing distribution
#'
#' @title Computing non-parametric mixing distribution
#' @param x a object from implemented family
#' @param mix initial mixing distribution (proper)
#' @param tol tolerance to stop the code
#' @param maxiter maximum iteration to stop the code
#' @param verbose logical. If TRUE, the intermediate results will be printed.
#' @export
computemixdist = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  UseMethod("computemixdist")
}

#' computing non-parametric mixing distribution with estimated proportion at 0
#'
#' This is a S3 generic function for computing non-parametric mixing
#' distribution with estimated proportion at 0. Different families will
#' have different threshold value.
#'
#' @title Computing non-parametric mixing distribution with estimated proportion at 0
#' @param x a object from implemented family
#' @param val the threshold value.
#' @param mix initial mixing distribution (proper)
#' @param tol tolerance to stop the code
#' @param maxiter maximum iteration to stop the code
#' @param verbose logical. If TRUE, the intermediate results will be printed.
#' @export
estpi0 = function(x, val, mix, tol, maxiter, verbose){
  UseMethod("estpi0")
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

checklossfunction = function(x, mu0, pi0, eta, p, maxit = 100){
  llorigin = lossfunction(x, mu0, pi0)
  con = sum(p * eta)
  sigma = 0.5; alpha = 1/3; u = -1;
  repeat{
    u = u + 1
    pi0new = pi0 + sigma^u * eta
    lhs = -lossfunction(x, mu0, pi0new)
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
}

collapsemix = function(x, mu0, pi0, tol = 1e-6){
  ll = lossfunction(x, mu0, pi0)
  ntol = max(tol * 0.1, ll * 1e-16)
  repeat{
    if (length(mu0) <= 1)
      break
    prec = c(10 * min(diff(mu0)))
    r = unique.disc(mu0, pi0, c(prec, 0))
    nll = lossfunction(x, r$pt, r$pr)
    if (nll <= ll + ntol){
      mu0 = r$pt; pi0 = r$pr
    }else{
      break
    }
  }
  list(pt = mu0, pr = pi0)
}
