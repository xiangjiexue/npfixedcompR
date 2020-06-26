#' @rdname makeobject
#' @export
makeobject.npnormllw = function(v, mu0, pi0, beta, order = -4){
  if (class(v) == "npnormllw"){
    # update information
    x = v
    if (!missing(mu0)) x$mu0 = mu0
    if (!missing(pi0)) x$pi0 = pi0
    if (!missing(beta)) x$beta = beta
    x$precompute = dnpdiscnorm(x$v, mu0 = x$mu0, pi0 = x$pi0, sd = x$beta, h = x$h)
  }

  if (is.numeric(v)){
    if (missing(mu0)) mu0 = 0
    if (missing(pi0)) pi0 = 0
    if (missing(beta)) beta = 1
    bindata = bin(v, order = order)
    x = list(v = bindata$v, w = bindata$w, mu0 = mu0, pi0 = pi0, beta = beta, h = 10^order,
             precompute = dnpdiscnorm(bindata$v, mu0 = mu0, pi0 = pi0, sd = beta, h = 10^order))
    attr(x, "class") = "npnormllw"
  }

  x
}

lossfunction.npnormllw = function(x, mu0, pi0){
  -sum(log(dnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h) + x$precompute) * x$w)
}

gradientfunction.npnormllw = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  if (!is.null(x$flexden)){
    flexden = x$flexden
  }else{
    flexden = dnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h)
  }
  fullden = flexden + x$precompute
  murep = rep(mu, rep(length(x$v), length(mu)))
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (order[1] == 1){
    temp = ddiscnorm(x$v, mean = murep, sd = x$beta, h = x$h) * sum(pi0)
    ans$d0 = .colSums((flexden - temp) / fullden * x$w, m = length(x$v), n = length(mu))
  }
  if (order[2] == 1){
    temp = dnorm(x$v + x$h, murep, sd = x$beta) - dnorm(x$v, murep, sd = x$beta)
    ans$d1 = .colSums(temp / fullden * x$w, m = length(x$v), n = length(mu)) * sum(pi0)
  }
  if (order[3] == 1){
    temp = dnorm(x$v + x$h, murep, sd = x$beta) * (x$v - murep) -
      dnorm(x$v, murep, sd = x$beta) * (x$v + x$h - murep)
    ans$d2 = .colSums(temp / fullden * x$w, m = length(x$v), n = length(mu)) / x$beta^2
  }

  ans
}

#' @export
computemixdist.npnormllw = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    breaks = pmax(ceiling(diff(rx) / (5 * x$beta)), 10)   # number of breaks
    r = whist(x$v, x$w, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    r$density = pmax(0, r$density  / sum(r$density) - pnpnorm(r$breaks[-1], mu0 = x$mu0, pi0 = x$pi0, sd = x$beta) +
                       pnpnorm(r$breaks[-length(r$breaks)], mu0 = x$mu0, pi0 = x$pi0, sd = x$beta))
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0] / sum(r$density)
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  x$fn = function(mu0, pi0) dnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h)
  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = gridpointsnpnorm(x)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol, method = "d1"))
    pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
    sp = ddiscnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$beta, h = x$h)
    fp = .rowSums(sp[1:(length(x$v) * length(mu0))] * rep(pi0, rep(length(x$v), length(pi0))), m = length(x$v), n = length(pi0)) + x$precompute
    S = sp / fp
    dim(S) = c(length(x$v), length(pi0new))
    a = 2 - x$precompute / fp
    newweight = pnnls(S * sqrt(x$w), a * sqrt(x$w), sum = 1 - sum(x$pi0))$x
    r = checklossfunction(x, mu0new, pi0new, newweight - pi0new, crossprod(S, x$w))
    r = collapsemix(x, r$pt, r$pr, tol)
    mu0 = r$pt; pi0 = r$pr
    iter = iter + 1
    nloss = lossfunction(x, mu0, pi0)

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

  mu0new = c(mu0, x$mu0); pi0new = c(pi0, x$pi0)
  index = order(mu0new, decreasing = FALSE)
  r = unique.disc(mu0new[index], pi0new[index])

  ans = list(iter = iter,
             family = "npnorm",
             max.gradient = min(gradientfunction(x, mu0, mu0, pi0, order = c(1, 0, 0))$d0),
             mix = list(pt = r$pt, pr = r$pr),
             ll = nloss + sum(x$w) * log(x$h),
             beta = x$beta,
             dd0 = gradientfunction(x, 0, mu0, pi0, order = c(1, 0, 0))$d0,
             h = x$h,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

#' @rdname estpi0
#' @export
estpi0.npnormllw = function(x, val = 0.5 * log(sum(x$w)), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x = makeobject(x, pi0 = 1 - tol / 2, method = attr(x, "class"))
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x = makeobject(x, pi0 = 0, method = attr(x, "class"))
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll - r0$ll < val){
    r = list(iter = 0,
             family = "npnorm",
             max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             ll = lossfunction(x, mu0 = 0, pi0 = 1) + sum(x$w) * log(x$h),
             dd0 = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             beta = x$beta,
             h = x$h,
             convergence = 0)
  }else{
    r = solveestpi0(x = x, init = dnpdiscnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$beta, h = x$h) / ddiscnorm(0, sd = x$beta, h = x$h),
                    val = -r0$ll - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
