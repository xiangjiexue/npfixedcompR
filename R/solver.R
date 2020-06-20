# In this case can f be a function with solving mu of order = c(0, 1, 0)
solvegradientsingle = function(x, mu0, pi0, lower = min(x$v), upper = max(x$v), tol = 1e-6){
  uniroot(function(mu_){
    gradientfunction(x = x, mu = mu_, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  }, lower = lower, upper = upper, tol = tol)$root
}

solvegradientmultiple = function(x, mu0, pi0, points, tol = 1e-6){
  pointsval = gradientfunction(x = x, mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
  if (length(index) >= 1){
    r = sapply(index, function(ddd){
      solvegradientsingle(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
    })
    r = c(min(x$v), r, max(x$v))
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = range(x$v)[gradientfunction(x = x, mu = range(x$v), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }
  r
}

# find the proportion of zero in the context of minimisation problem
# pi0 = 0 should always be negative
solveestpi0 = function(x, init, val, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  neg = 0; pos = 1
  iter = 1;
  x = makeobject(x, mu0 = x$mu0, pi0 = init, beta = x$sd, method = attr(x, "class"))
  r = computemixdist(x, mix = mix, maxiter = maxiter, tol = tol)
  d1 = r$ll + val
  d2 = r$dd0
  repeat{
    if (verbose){
      cat("Iteration", iter, "\n")
      cat("Current pi0 value =", round(init, -ceiling(log10(tol))), "with function value", d1,
          "(iterations", r$iter, ")\n")
      cat("Current Bracket: lower =", neg, "upper =", pos, "\n")
    }

    if (iter > maxiter | abs(d1) < tol | pos - neg < tol)
      break
    iter = iter + 1
    if (d1 < 0){
      neg = max(neg, init)
    }else{
      pos = min(pos, init)
    }

    init = init - d1 / d2

    if ((init - neg) < 0.1 * (pos - neg) | init < neg){
      init = (pos - neg) * 0.1 + neg
    }

    if ((init - neg) > 0.9 * (pos - neg) | init > pos){
      init = (pos - neg) * 0.9 + neg
    }

    x = makeobject(x, mu0 = x$mu0, pi0 = init, beta = x$sd, method = attr(x, "class"))
    j0 = r$mix$pt == 0
    newmix = list(pt = r$mix$pt[!j0], pr = r$mix$pr[!j0] / sum(r$mix$pr[!j0]))
    r = computemixdist(x, mix = newmix, maxiter = maxiter, tol = tol)
    d1 = r$ll + val
    d2 = r$dd0
  }

  r
}
