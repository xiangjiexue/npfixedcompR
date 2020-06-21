#' @rdname makeobject
#' @export
makeobject.npnormcvmw = function(v, mu0 = 0, pi0 = 0, beta = 1, order = -4){
  if (class(v) == "npnormcvmw"){
    # update information
    x = v
    x$mu0 = mu0; x$pi0 = pi0; x$sd = beta
    x$precompute = cumsum(x$w) / sum(x$w) - 0.5 / sum(x$w) -
      pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = beta, h = x$h)
  }

  if (is.numeric(v)){
    bindata = bin(v, order = order)
    x = list(v = bindata$v, w = bindata$w, mu0 = mu0, pi0 = pi0, sd = beta, h = 10^order,
             precompute = cumsum(bindata$w) / sum(bindata$w) - 0.5 / sum(bindata$w) -
               pnpdiscnorm(bindata$v, mu0 = mu0, pi0 = pi0, sd = beta, h = 10^order))
    attr(x, "class") = "npnormcvmw"
  }

  x
}

lossfunction.npnormcvmw = function(x, mu0, pi0){
  sum((pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$sd, h = x$h) - x$precompute)^2 * x$w)
}

gradientfunction.npnormcvmw = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  flexden = pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$sd, h = x$h)
  fullden = flexden - x$precompute
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (order[1] == 1){
    temp = pdiscnorm(x$v, mean = rep(mu, rep(length(x$v), length(mu))), sd = x$sd, h = x$h) * sum(pi0)
    ans$d0 = .colSums((temp - flexden) * fullden * x$w, m = length(x$v), n = length(mu))
  }
  if (any(order[2:3] == 1)){
    temp = dnorm(x$v + x$h, rep(mu, rep(length(x$v), length(mu))), sd = x$sd) * fullden * x$w
  }
  if (order[2] == 1){
    ans$d1 = .colSums(temp, m = length(x$v), n = length(mu)) * -2 * sum(pi0)
  }
  if (order[3] == 1){
    xminusmu = x$v + x$h - rep(mu, rep(length(x$v), length(mu)))
    ans$d2 = .colSums(temp * xminusmu, m = length(x$v), n = length(mu)) * -2 * sum(pi0)
  }

  ans
}

#' @export
computemixdist.npnormcvmw = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    breaks = pmax(ceiling(diff(rx) / (5 * x$sd)), 5)   # number of breaks
    r = whist(x$v, x$w, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    r$density = pmax(0, r$density - pnpnorm(r$breaks[-1], mu0 = x$mu0, pi0 = x$pi0, sd = x$sd) +
                       pnpnorm(r$breaks[-length(r$breaks)], mu0 = x$mu0, pi0 = x$pi0, sd = x$sd))
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0] / sum(r$density)
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = c(min(x$v) - 3 * x$sd, gridpointsnpnorm(x), max(x$v) + 3 * x$sd)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol))
    S = pdiscnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$sd, h = x$h)
    dim(S) = c(length(x$v), length(mu0new))
    pi0new = pnnls(S * sqrt(x$w), x$precompute * sqrt(x$w), sum = 1 - sum(x$pi0))$x
    r = collapsemix(x, mu0new, pi0new, tol)
    mu0 = r$pt; pi0 = r$pr
    iter = iter + 1
    nloss = lossfunction(x, mu0, pi0)

    if (verbose){
      cat("Iteration: ", iter, "\n")
      cat(paste0("Support Point ", round(r$pt, -ceiling(log10(tol))), " with probability ", round(r$pr, -ceiling(log10(tol))), "\n"))
      cat("Current Cramer-von Mises Loss ", as.character(round(nloss, -ceiling(log10(tol)))), "\n")
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
             ll = nloss,
             dd0 = gradientfunction(x, 0, r$pt, r$pr, order = c(1, 0, 0))$d0,
             h = x$h,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

#' @rdname estpi0
#' @export
estpi0.npnormcvmw = function(x, val = qCvM(0.05, lower.tail = FALSE), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x = makeobject(x, mu0 = 0, pi0 = 1 - tol / 2, beta = x$sd, method = attr(x, "class"))
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x = makeobject(x, mu0 = 0, pi0 = 0, beta = x$sd, method = attr(x, "class"))
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  nval = sum(x$w) / 3 - sum(x$w * (cumsum(x$w) / sum(x$w) - 0.5 / sum(x$w))^2)

  if (r1$ll + nval < val){
    r = list(iter = 0,
             family = "npnorm",
             max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             ll = lossfunction(x, mu0 = 0, pi0 = 1),
             dd0 = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             convergence = 0)
  }else{
    r = solveestpi0(x = x, init = dnpdiscnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$sd, h = x$h) / ddiscnorm(0, sd = x$sd, h = x$h),
                    val = nval - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
