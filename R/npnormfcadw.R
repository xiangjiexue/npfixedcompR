#' @rdname makeobject
#' @export
makeobject.npnormadw = function(v, mu0, pi0, beta, order = -4){
  if (class(v) == "npnormadw"){
    # update information
    x = v
    if (!missing(mu0)) x$mu0 = mu0
    if (!missing(pi0)) x$pi0 = pi0
    if (!missing(beta)) x$beta = beta
    x$precompute1 = pnpdiscnorm(v$v, mu0 = x$mu0, pi0 = x$pi0, sd = x$beta, h = x$h)
    x$precompute2 = pnpdiscnorm(v$v, mu0 = x$mu0, pi0 = x$pi0, sd = x$beta, lower.tail = FALSE, h = x$h)
  }

  if (is.numeric(v)){
    bindata = bin(v, order = order)
    if (missing(mu0)) mu0 = 0
    if (missing(pi0)) pi0 = 0
    if (missing(beta)) beta = 1
    x = list(v = bindata$v, w = bindata$w, mu0 = mu0, pi0 = pi0, beta = beta, h = 10^order,
             a1 = (2 * cumsum(bindata$w) - 1) / sum(bindata$w) * bindata$w,
             a2 = (2 + (-2 * cumsum(bindata$w) + 1) / sum(bindata$w)) * bindata$w,
             precompute1 = pnpdiscnorm(bindata$v, mu0 = mu0, pi0 = pi0, sd = beta, h = 10^order),
             precompute2 = pnpdiscnorm(bindata$v, mu0 = mu0, pi0 = pi0, sd = beta, lower.tail = FALSE, h = 10^order))
    attr(x, "class") = "npnormadw"
  }

  x
}

lossfunction.npnormadw = function(x, mu0, pi0){
  temp = pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h) + x$precompute1
  -sum(log(temp) * x$a1 + log1p(-temp) * x$a2)
}

gradientfunction.npnormadw = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  if (!is.null(x$flexden)){
    flexden = x$flexden$p1; flexden2 = x$flexden$p2
  }else{
    flexden = pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h)
    flexden2 = pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, lower.tail = FALSE, h = x$h)
  }
  murep = rep(mu, rep(length(x$v), length(mu)))
  fullden = flexden + x$precompute1
  fullden2 = flexden2 + x$precompute2
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (order[1] == 1){
    addconst = -sum(x$precompute1 * x$a1 / fullden + x$precompute2 * x$a2 / fullden2) + 2 * sum(x$w)
    ans$d0 = .colSums(pdiscnorm(x$v, murep, sd = x$beta, h = x$h) * x$a1 / fullden +
                        pdiscnorm(x$v, murep, sd = x$beta, lower.tail = FALSE, h = x$h) * x$a2 / fullden2,
                      m = length(x$v), n = length(mu)) * -sum(pi0) + addconst
  }
  if (any(order[2:3] == 1)){
    temp = dnorm(x$v + x$h, mean = murep, sd = x$beta) * (x$a1 / fullden - x$a2 / fullden2)
  }
  if (order[2] == 1){
    ans$d1 = .colSums(temp, m = length(x$v), n = length(mu)) * sum(pi0)
  }
  if (order[3] == 1){
    ans$d2 = .colSums((x$v + x$h - murep) * temp, m = length(x$v), n = length(mu)) * sum(pi0) / x$beta^2
  }

  ans
}

#' @export
computemixdist.npnormadw = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
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

  x$fn = function(mu0, pi0) list(p1 = pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h),
                                 p2 = pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h, lower.tail = FALSE))

  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = c(min(x$v) - 3 * x$beta, gridpointsnpnorm(x), max(x$v) + 3 * x$beta)
  con1 = x$precompute1 * x$a1; con2 = x$precompute2 * x$a2

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol, method = "d2"))
    pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
    sf = pdiscnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$beta, h = x$h)
    ss = .rowSums(sf[1:(length(x$v) * length(mu0))] * rep(pi0, rep(length(x$v), length(pi0))), m = length(x$v), n = length(pi0)) + x$precompute1
    S = sf / ss
    dim(S) = c(length(x$v), length(pi0new))

    uf = pdiscnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$beta, lower.tail = FALSE, h = x$h)
    us = .rowSums(uf[1:(length(x$v) * length(mu0))] * rep(pi0, rep(length(x$v), length(pi0))), m = length(x$v), n = length(pi0)) + x$precompute2
    U = uf / us
    dim(U) = c(length(x$v), length(pi0new))

    S1 = crossprod(S, S * x$a1) + crossprod(U, U * x$a2)
    S2 = colSums(S * x$a1) + colSums(U * x$a2)

    newweight = pnnqp(S1, -2 * S2 + crossprod(S, con1 / ss) + crossprod(U, con2 / us), sum = 1 - sum(x$pi0))$x
    r = checklossfunction(x, mu0new, pi0new, newweight - pi0new, S2)
    r = collapsemix(x, r$pt, r$pr, tol)
    mu0 = r$pt; pi0 = r$pr
    iter = iter + 1
    nloss = lossfunction(x, mu0, pi0)

    if (verbose){
      cat("Iteration: ", iter, "\n")
      cat(paste0("Support Point ", round(r$pt, -ceiling(log10(tol))), " with probability ", round(r$pr, -ceiling(log10(tol))), "\n"))
      cat("Current Anderson-Darling Loss ", as.character(round(-nloss, -ceiling(log10(tol)))), "\n")
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
             beta = x$beta,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

estpi0d.npnormadw = function(x, mu0, pi0){
  ans = vector("list", 2)
  names(ans) = c("d2", "d3")
  S1 = pdiscnorm(x$v, sd = x$beta, h = x$h) / pnpdiscnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta, h = x$h) - 1
  S2 = pdiscnorm(x$v, sd = x$beta, lower.tail = FALSE, h = x$h) / pnpdiscnorm(x$v, mu0, pi0, sd = x$beta, lower.tail = FALSE, h = x$h) - 1
  ans$d2 = -sum((S1 * x$a1 + S2 * x$a2) * x$w); ans$d3 = sum((S1^2 * x$a1 + S2^2 * x$a2) * x$w)
  ans
}

#' @rdname estpi0
#' @export
estpi0.npnormadw = function(x, val = qAD(0.05, lower.tail = FALSE), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x = makeobject(x, pi0 = 1 - tol / 2, method = attr(x, "class"))
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x = makeobject(x, pi0 = 0, method = attr(x, "class"))
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll - sum(x$w) < val){
    r = list(iter = 0,
             family = "npnorm",
             max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             ll = lossfunction(x, mu0 = 0, pi0 = 1),
             beta = x$beta,
             convergence = 0)
  }else{
    r = solveestpi0(x = x, init = dnpdiscnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$beta, h = x$h) / ddiscnorm(0, sd = x$beta, h = x$h),
                    val = -sum(x$w) - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
