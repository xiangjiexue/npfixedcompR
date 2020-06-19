#' npfixedcompR
#'
#' npfixedcompR
#'
#' @docType package
#' @author Xiangjie Xue
#' @importFrom stats dnorm uniroot
#' @importFrom lsei pnnls
#' @name npfixedcompR
#' NULL

lossfunction = function(x, mu0, pi0){
  UseMethod("lossfunction", x)
}


#' gradient function
#'
#' This is a S3 generic function returning gradident function with order specified as a list
gradientfunction = function(x, mu, mu0, pi0, order){
  UseMethod("gradientfunction", x)
}

computemixdist = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  UseMethod("computemixdist")
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
