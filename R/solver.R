# In this case can f be a function with solving mu of order = c(0, 1, 0)
solvegradientsingled1 = function(x, mu0, pi0, lower = min(x$v), upper = max(x$v), tol = 1e-6){
  uniroot(function(mu_){
    gradientfunction(x = x, mu = mu_, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  }, lower = lower, upper = upper, tol = tol)$root
}

solvegradientmultipled1 = function(x, mu0, pi0, points, tol = 1e-6){
  pointsval = gradientfunction(x = x, mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
  if (length(index) >= 1){
    r = sapply(index, function(ddd){
      solvegradientsingled1(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
    })
    r = c(min(x$v), r, max(x$v))
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = range(x$v)[gradientfunction(x = x, mu = range(x$v), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }
  r
}

solvegradientsingled0 = function(x, mu0, pi0, lower = min(x$v), upper = max(x$v), tol = 1e-6){
  optimise(function(mu_){
    gradientfunction(x = x, mu = mu_, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0
  }, lower = lower, upper = upper, tol = tol)$minimum
}

solvegradientmultipled0 = function(x, mu0, pi0, points, tol = 1e-6){
  pointsval = diff(gradientfunction(x = x, mu = points, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0)
  index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
  if (length(index) >= 1){
    r = sapply(index, function(ddd){
      solvegradientsingled0(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
    })
    r = c(min(x$v), r, max(x$v))
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = range(x$v)[gradientfunction(x = x, mu = range(x$v), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }
  r
}

solvegradientsingled2 = function(x, mu0, pi0, lower = min(x$v), upper = max(x$v), tol = 1e-6){
  neg = lower; pos = upper; init = mean(c(neg, pos))
  negval = NA; posval = NA
  r = gradientfunction(x, init, mu0, pi0, order = c(0, 1, 1))
  d1 = r$d1
  d2 = r$d2
  repeat{
    if (abs(d1) < tol | pos - neg < tol)
      break
    if (d1 < 0){
      neg = max(neg, init)
      negval = max(negval, d1, na.rm = TRUE)
    }else{
      pos = min(pos, init)
      posval = min(posval, d1, na.rm = TRUE)
    }

    init = init - d1 / d2

    if ((init - neg) < 0.1 * (pos - neg) | init < neg){
      if (!is.na(negval) & !is.na(posval)){
        init = (1 + negval / (posval - negval)) * neg + (1 - posval / (posval - negval)) * pos
      }else{
        init = (pos - neg) * 0.1 + neg
      }

    }

    if ((init - neg) > 0.9 * (pos - neg) | init > pos){
      if (!is.na(negval) & !is.na(posval)){
        init = (1 + negval / (posval - negval)) * neg + (1 - posval / (posval - negval)) * pos
      }else{
        init = (pos - neg) * 0.9 + neg
      }
    }

    r = gradientfunction(x, init, mu0, pi0, order = c(0, 1, 1))
    d1 = r$d1
    d2 = r$d2
  }

  init
}

solvegradientmultipled2 = function(x, mu0, pi0, points, tol = 1e-6){
  pointsval = gradientfunction(x = x, mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
  if (length(index) >= 1){
    r = sapply(index, function(ddd){
      solvegradientsingled2(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
    })
    r = c(min(x$v), r, max(x$v))
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = range(x$v)[gradientfunction(x = x, mu = range(x$v), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }
  r
}

solvegradientmultiple = function(x, mu0, pi0, points, tol = 1e-6, method = "auto"){
  switch(method,
         "auto" = {
          # minimal sample to test implementation
          xx = makeobject(1, method = attr(x, "class"))
          rr = gradientfunction(xx, 1, mu0 = 0, pi0 = 1, order = c(1, 1, 1))
          if (!is.null(rr$d2)){
            solvegradientmultipled2(x, mu0, pi0, points, tol)
          }else{
            if (!is.null(rr$d1)){
              solvegradientmultipled1(x, mu0, pi0, points, tol)
            }else{
              solvegradientmultipled0(x, mu0, pi0, points, tol)
            }
          }
         },
         "d1" = solvegradientmultipled1(x, mu0, pi0, points, tol),
         "d2" = solvegradientmultipled2(x, mu0, pi0, points, tol),
         "d0" = solvegradientmultipled0(x, mu0, pi0, points, tol))
}

# find the proportion of zero in the context of minimisation problem
# pi0 = 0 should always be negative
solveestpi0 = function(x, init, val, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  neg = 0; pos = 1 - tol / 2
  negval = NA; posval = NA
  iter = 1;
  x = makeobject(x, pi0 = init, method = attr(x, "class"))
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
      negval = max(negval, d1, na.rm = TRUE)
    }else{
      pos = min(pos, init)
      posval = min(posval, d1, na.rm = TRUE)
    }

    init = init - d1 / d2

    if ((init - neg) < 0.1 * (pos - neg) | init < neg){
      if (!is.na(negval) & !is.na(posval)){
        init = (1 + negval / (posval - negval)) * neg + (1 - posval / (posval - negval)) * pos
      }else{
        init = (pos - neg) * 0.1 + neg
      }

    }

    if ((init - neg) > 0.9 * (pos - neg) | init > pos){
      if (!is.na(negval) & !is.na(posval)){
        init = (1 + negval / (posval - negval)) * neg + (1 - posval / (posval - negval)) * pos
      }else{
        init = (pos - neg) * 0.9 + neg
      }
    }

    x = makeobject(x, pi0 = init, method = attr(x, "class"))
    j0 = r$mix$pt == 0
    newmix = list(pt = r$mix$pt[!j0], pr = r$mix$pr[!j0] / sum(r$mix$pr[!j0]))
    r = computemixdist(x, mix = newmix, maxiter = maxiter, tol = tol)
    d1 = r$ll + val
    d2 = r$dd0
  }

  r
}
