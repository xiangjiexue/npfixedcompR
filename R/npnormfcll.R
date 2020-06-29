#' @rdname makeobject
#' @export
makeobject.npnormll = function(v, mu0, pi0, beta){
  if (class(v) == "npnormll"){
    x = v
    if (!missing(mu0)) x$mu0 = mu0
    if (!missing(pi0)) x$pi0 = pi0
    if (!missing(beta)) x$beta = beta
    # update information
    x$precompute = dnpnorm(x$v, mu0 = x$mu0, pi0 = x$pi0, sd = x$beta)
  }

  if (is.numeric(v)){
    if (missing(mu0)) mu0 = 0
    if (missing(pi0)) pi0 = 0
    if (missing(beta)) beta = 1
    x = list(v = v, mu0 = mu0, pi0 = pi0, beta = beta,
             precompute = dnpnorm(v, mu0 = mu0, pi0 = pi0, sd = beta))
    attr(x, "class") = "npnormll"
  }

  x
}

lossfunction.npnormll = function(x, mu0, pi0){
  -sum(log(dnpnorm(x$v, mu0 = mu0, pi0 = pi0, x$beta) + x$precompute))
}

gradientfunction.npnormll = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  if (!is.null(x$flexden)){
    flexden = x$flexden
  }else{
    flexden = dnpnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta)
  }
  murep = rep(mu, rep(length(x$v), length(mu)))
  temp = dnorm(x$v, mean = murep, sd = x$beta) * sum(pi0)
  fullden = flexden + x$precompute
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (order[1] == 1){
    ans$d0 = .colSums((flexden - temp) / fullden, m = length(x$v), n = length(mu))
  }
  if (any(order[2:3] == 1)){
    xminusmu = x$v - murep
    temp2 = temp / fullden
  }
  if (order[2] == 1){
    ans$d1 = .colSums(temp2 * xminusmu, m = length(x$v), n = length(mu)) / -x$beta^2
  }
  if (order[3] == 1){
    ans$d2 = .colSums(temp2 * (xminusmu^2 - x$beta^2), m = length(x$v), n = length(mu)) / x$beta^4
  }

  ans
}

#' @export
computemixdist.npnormll = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    breaks = pmax(ceiling(diff(rx) / (5 * x$beta)), 10)   # number of breaks
    r = whist(x$v, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    r$density = pmax(0, r$density  / sum(r$density) - pnpnorm(r$breaks[-1], mu0 = x$mu0, pi0 = x$pi0, sd = x$beta) +
                       pnpnorm(r$breaks[-length(r$breaks)], mu0 = x$mu0, pi0 = x$pi0, sd = x$beta))
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0] / sum(r$density)
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  x$fn = function(mu0, pi0) dnpnorm(x$v, mu0, pi0, x$beta)

  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = gridpointsnpnorm(x)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol, method = "d1"))
    pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
    sp = dnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$beta)
    fp = .rowSums(sp[1:(length(x$v) * length(mu0))] * rep(pi0, rep(length(x$v), length(pi0))), m = length(x$v), n = length(pi0)) + x$precompute
    S = sp / fp
    dim(S) = c(length(x$v), length(pi0new))
    a = 2 - x$precompute / fp
    newweight = pnnls(S, a, sum = 1 - sum(x$pi0))$x
    r = checklossfunction(x, mu0new, pi0new, newweight - pi0new, colSums(S))
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
             ll = nloss,
             beta = x$beta,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

estpi0d.npnormll = function(x, mu0, pi0){
  ans = vector("list", 2)
  names(ans) = c("d2", "d3")
  S = dnorm(x$v, sd = x$beta) / dnpnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$beta) - 1
  ans$d2 = -sum(S); ans$d3 = sum(S^2)
  ans
}

#' @rdname estpi0
#' @export
estpi0.npnormll = function(x, val = 0.5 * log(length(x$v)), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x = makeobject(x, pi0 = 1 - tol / 2, method = attr(x, "class"))
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x = makeobject(x, pi0 = 0, method = attr(x, "class"))
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll - r0$ll < val){
    r = list(iter = 0,
             family = "npnorm",
             max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             beta = x$beta,
             ll = lossfunction(x, mu0 = 0, pi0 = 1),
             convergence = 0)
  }else{
    r = solveestpi0(x = x, init = dnpnorm(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, sd = x$beta) * sqrt(2 * base::pi) * x$beta,
                    val = -r0$ll - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
