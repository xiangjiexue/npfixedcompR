lossfunction.npnormfcll = function(x, mu0, pi0){
  -sum(log(dnpnorm(x$v, mu0 = mu0, pi0 = pi0, x$sd) + x$precompute))
}

gradientfunction.npnormfcll = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  flexden = dnpnorm(x$v, mu0 = mu0, pi0 = pi0, sd = x$sd)
  temp = dnorm(x$v, mean = rep(mu, rep(length(x$v), length(mu))), sd = x$sd) * sum(pi0)
  fullden = flexden + x$precompute
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (any(order[2:3] == 1)){
    xminusmu = rep(mu, rep(length(x$v), length(mu))) - x$v
  }
  if (order[1] == 1){
    ans$d0 = .colSums((flexden - temp) / fullden, m = length(x$v), n = length(mu))
  }
  if (order[2] == 1){
    ans$d1 = .colSums(temp * xminusmu / fullden, m = length(x$v), n = length(mu)) / x$sd^2
  }
  if (order[3] == 1){
    ans$d2 = .colSums(temp * (xminusmu^2 - x$sd^2) / fullden, m = length(x$v), n = length(mu)) / x$sd^4
  }

  ans
}

computemixdist.npnormfcll = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    breaks = pmax(ceiling(diff(rx) / (5*x$sd)), 5)   # number of breaks
    r = hist(x$v, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0]
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  pi0 = pi0 * (1 - sum(x$pi0))
  x$precompute = dnpnorm(x$v, mu0 = x$mu0, pi0 = x$pi0, sd = x$sd)

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = gridpointsnpnorm(x)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol))
    pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
    sp = dnorm(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), sd = x$sd)
    fp = .rowSums(sp[1:(length(x$v) * length(mu0))] * rep(pi0, rep(length(x$v), length(pi0))), m = length(x$v), n = length(pi0)) + x$precompute
    S = sp / fp
    dim(S) = c(length(x$v), length(pi0new))
    a = 2 - x$precompute / fp
    newweight = pnnls(S, a, sum = 1 - sum(x$pi0))$x
    r = checklossfunction(x, mu0new, pi0new, newweight - pi0new, colSums(S) / (1 - sum(x$pi0)))
    r = collapsemix(x, r$pt, r$pr, tol)
    mu0 = r$pt; pi0 = r$pr
    iter = iter + 1
    nloss = lossfunction(x, mu0, pi0)

    if (verbose){
      cat(paste("Support Point", round(r$pt, 3), "with probability", round(r$pr, 3), "\n"))
      cat("Current log-likelihood", -nloss, "\n")
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
             max.gradient = -min(gradientfunction(x, mu0, mu0, pi0, order = c(1, 0, 0))$d0),
             mix = list(pt = r$pt, pr = r$pr),
             ll = -nloss,
             dd0 = -gradientfunction(x, 0, r$pt, r$pr, order = c(1, 0, 0))$d0,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}
