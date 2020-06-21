#' @rdname makeobject
#' @export
makeobject.npnormcvm = function(v, mu0, pi0, beta){
  if (class(v) == "npnormcvm"){
    # update information
    x = v
    if (!missing(mu0)) x$mu0 = mu0
    if (!missing(pi0)) x$pi0 = pi0
    if (!missing(beta)) x$beta = beta
    x$precompute = seq(from = 0.5 / length(x$v), to = 1 - 0.5 / length(x$v), length = length(x$v)) -
      pnpnorm(x$v, mu0 = mu0, pi0 = pi0, sd = beta)
  }

  if (is.numeric(v)){
    v = sort(v, decreasing = FALSE)
    if (missing(mu0)) mu0 = 0
    if (missing(pi0)) pi0 = 0
    if (missing(beta)) beta = 1
    x = list(v = v, mu0 = mu0, pi0 = pi0, beta = beta,
             precompute = seq(from = 0.5 / length(v), to = 1 - 0.5 / length(v), length = length(v)) -
               pnpnorm(v, mu0 = mu0, pi0 = pi0, sd = beta))
    attr(x, "class") = "npnormcvm"
  }

  x
}

lossfunction.npnormcvm = function(x, mu0, pi0){
  sum((pnpnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta) - x$precompute)^2)
}

gradientfunction.npnormcvm = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  flexden = pnpnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta)
  fullden = flexden - x$precompute
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (order[1] == 1){
    temp = pnorm(x$v, mean = rep(mu, rep(length(x$v), length(mu))), sd = x$beta) * sum(pi0)
    ans$d0 = .colSums((temp - flexden) * fullden, m = length(x$v), n = length(mu))
  }
  if (any(order[2:3] == 1)){
    temp = dnorm(x$v, mean = rep(mu, rep(length(x$v), length(mu))), sd = x$beta) * sum(pi0)
  }
  if (order[2] == 1){
    ans$d1 = .colSums(temp * fullden, m = length(x$v), n = length(mu)) * -2
  }
  if (order[3] == 1){
    xminusmu = x$v - rep(mu, rep(length(x$v), length(mu)))
    ans$d2 = .colSums(temp * xminusmu * fullden, m = length(x$v), n = length(mu)) * -2
  }

  ans
}

#' @export
computemixdist.npnormcvm = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    breaks = pmax(ceiling(diff(rx) / (5 * x$beta)), 5)   # number of breaks
    r = whist(x$v, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    r$density = pmax(0, r$density - pnpnorm(r$breaks[-1], mu0 = x$mu0, pi0 = x$pi0, sd = x$beta) +
                       pnpnorm(r$breaks[-length(r$breaks)], mu0 = x$mu0, pi0 = x$pi0, sd = x$beta))
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0] / sum(r$density)
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = c(min(x$v) - 3 * x$beta, gridpointsnpnorm(x), max(x$v) + 3 * x$beta)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol, method = "d2"))
    S = pnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$beta)
    dim(S) = c(length(x$v), length(mu0new))
    pi0new = pnnls(S, x$precompute, sum = 1 - sum(x$pi0))$x
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
             beta = x$beta,
             dd0 = gradientfunction(x, 0, mu0, pi0, order = c(1, 0, 0))$d0,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

#' @rdname estpi0
#' @export
estpi0.npnormcvm = function(x, val = qCvM(0.05, lower.tail = FALSE), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x = makeobject(x, pi0 = 1 - tol / 2, method = attr(x, "class"))
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x = makeobject(x, pi0 = 0, method = attr(x, "class"))
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll + 1/ 12 / length(x$v) < val){
    r = list(iter = 0,
             family = "npnorm",
             max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             ll = lossfunction(x, mu0 = 0, pi0 = 1),
             beta = x$beta,
             dd0 = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             convergence = 0)
  }else{
    r = solveestpi0(x = x, init = dnpnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$beta) * sqrt(2 * base::pi) * x$beta,
                    val = 1/ 12 / length(x$v) - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
