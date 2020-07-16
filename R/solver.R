# In this case can f be a function with solving mu of order = c(0, 1, 0)
solvegradientsingled1 = function(x, mu0, pi0, lower = min(x$v), upper = max(x$v), tol = 1e-6){
  # vectorised implementation
  a = lower; b = upper
  f1 = gradientfunction(x = x, mu = c(a, b), mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  fa = f1[1:length(a)]; fb = f1[(length(a) + 1) : length(f1)]

  fs = fb; s = b
  fc = numeric(length(a))
  index = rep(FALSE, length(a))
  repeat{
    if (sum(index) == length(a))
      break

    index = abs(fb) < tol | abs(fs) < tol | abs(b - a) < tol

    c1 = (a + b) / 2
    fc[!index] = gradientfunction(x, c1[!index], mu0, pi0, order = c(0, 1, 0))$d1

    # does not matter how the following changes for a[index], b[index], etc. changes, not referenced as long as
    # s is not assigned to a or b
    j0 = fa != fc & fb != fc
    s[j0] = a[j0] * fb[j0] * fc[j0] / (fa[j0] - fb[j0]) / (fa[j0] - fc[j0]) +
      b[j0] * fa[j0] * fc[j0] / (fb[j0] - fa[j0]) / (fb[j0] - fc[j0]) +
      c1[j0] * fa[j0] * fb[j0] / (fc[j0] - fa[j0]) / (fc[j0] - fb[j0])
    s[!j0] = b[!j0] - fb[!j0] * (b[!j0] - a[!j0]) / (fb[!j0] - fa[!j0])

    fs[!index] = gradientfunction(x, s[!index], mu0, pi0, order = c(0, 1, 0))$d1

    j0 = fs > 0
    j1 = s > a
    a[!j0 & !index & j1] = s[!j0 & !index & j1]
    fa[!j0 & !index & j1] = fs[!j0 & !index & j1]
    j1 = s < b
    b[j0 & !index & j1] = s[j0 & !index & j1]
    fb[j0 & !index & j1] = fs[j0 & !index & j1]

    j0 = fc > 0
    a[!j0 & !index] = c1[!j0 & !index]
    fa[!j0 & !index] = fc[!j0 & !index]
    b[j0 & !index] = c1[j0 & !index]
    fb[j0 & !index] = fc[j0 & !index]
  }

  b
}


solvegradientmultipled1 = function(x, mu0, pi0, points, tol = 1e-6){
  pointsval = gradientfunction(x = x, mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
  if (length(index) >= 1){
    # r = sapply(index, function(ddd){
    #   solvegradientsingled1(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
    # })
    r = c(points[1], solvegradientsingled1(x = x, mu0 = mu0, pi0 = pi0, lower = points[index], upper = points[index + 1], tol), points[length(points)])
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = c(points[1], points[length(points)])[gradientfunction(x = x, mu = c(points[1], points[length(points)]), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
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
    r = c(points[1], r, points[length(points)])
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = c(points[1], points[length(points)])[gradientfunction(x = x, mu = c(points[1], points[length(points)]), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }
  r
}

solvegradientsingled2 = function(x, mu0, pi0, lower = min(x$v), upper = max(x$v), tol = 1e-6){
  # This implements the vectorised version
  neg = lower; pos = upper; init = (neg + pos) / 2
  r = gradientfunction(x, init, mu0, pi0, order = c(0, 1, 1))
  d1 = r$d1
  d2 = r$d2
  d = numeric(length(lower))
  mflag = rep(FALSE, length(neg)) # There is no need to calculate the gradient function for mid-point again.

  index = rep(FALSE, length(neg))
  repeat{
    if (sum(index) == length(neg))
      break
    index = abs(d1) < tol | pos - neg < tol

    if (sum(mflag) > 0){
      s = (pos + neg) / 2
      d[!index & mflag] = gradientfunction(x, s[!index & mflag], mu0, pi0, order = c(0, 1, 0))$d1

      j0 = d < 0
      neg[j0 & !index & mflag] = s[j0 & !index & mflag]
      pos[!j0 & !index & mflag] = s[!j0 & !index & mflag]
    }

    j0 = d1 < 0
    neg[j0] = pmax(neg[j0], init[j0])
    pos[!j0] = pmin(pos[!j0], init[!j0])

    init[!index] = init[!index] - d1[!index] / d2[!index]
    mflag = init < neg | init > pos
    init[!index & mflag] = (pos[!index & mflag] + neg[!index & mflag]) / 2

    r = gradientfunction(x, init[!index], mu0, pi0, order = c(0, 1, 1))
    d1[!index] = r$d1
    d2[!index] = r$d2
  }

  init
}

solvegradientmultipled2 = function(x, mu0, pi0, points, tol = 1e-6){
  pointsval = gradientfunction(x = x, mu = points, mu0 = mu0, pi0 = pi0, order = c(0, 1, 0))$d1
  index = seq(1, by = 1, length = length(points) - 1)[pointsval[-length(pointsval)] < 0 & pointsval[-1] > 0]
  if (length(index) >= 1){
    r = c(points[1], solvegradientsingled2(x = x, mu0 = mu0, pi0 = pi0, lower = points[index], upper = points[index + 1], tol), points[length(points)])
    r = r[gradientfunction(x = x, mu = r, mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }else{
    r = c(points[1], points[length(points)])[gradientfunction(x = x, mu = c(points[1], points[length(points)]), mu0 = mu0, pi0 = pi0, order = c(1, 0, 0))$d0 < 0]
  }
  r
}

#' This function is used for finding the support points by finding the local minima.
#'
#' This series functions has three layers: The first layer \code{solvegradientmultiple} chooses the algorithm
#' to use, the second layer \code{solvegradientmultipled0}, \code{solvegradientmultipled1},
#' \code{solvegradientmultipled2} doing the simple derivative tests to narrow down the search space to be
#' the local minima (since all the optimisation problems have been done through minimisation), the third
#' layer \code{solvegradientsingled0}, \code{solvegradientsingled1}, \code{solvegradientsingled2} is the
#' larbour functions do the minimisation.
#'
#' Three different algorithms are implmented:
#'
#' - The derivative-free minimisation: the derivative-free minimisation is used when the gradient
#' function only have \code{d0} implemented. The implementation is the \code{\link{optimise}} function
#' in base R pacakge.
#'
#' - The derivative-free root-finding algorithm: the derivative-free root-finding is used when
#' the gradient function have \code{d1} implemented. The implementaion is the vectorised reduction version of the
#' modified Brent's method by Zhang (2011) with modifications tailored to this problem. The modification
#' ensure that in each iteration, the search space is at least half of the one from last iteration.
#'
#' - The first-order root-finding algorithm: the first-order root-finding is used when the gradient
#' function have both \code{d1} and \code{d2} implemented. The implementation is the vectorised reduction
#' Newton-Raphson method with modifications tailored to this problem. The modification ensures that
#' in each iteration, the search space is at least half of the one from last iteration.
#'
#' Vectorised reduction algorithms means that the function can have vectors of \code{lower} and \code{upper}
#' as inputs and when any element in each iteration satisfies the stopping criterion, it will not be
#' computed further. This should be more efficient when dealing with large datasets since logical operators
#' on a smaller set should always be computationally cheaper than keeping computing that support points to
#' an unnecessary accuracy.
#'
#' Note: This function is not an exported object.
#'
#' @title computing the support points
#' @param x a object from implemented family
#' @param mu0 a vector of support points
#' @param pi0 a vector of weights corresponding to the support points
#' @param points given grid points for finding the local minima (if exists)
#' @param tol tolerance
#' @param method the character strings specifying which algorithm to use. \code{method = "auto"} choose
#' the algorithm automatically, \code{method = "d0"} chooses the derivative-free minimisation,
#' \code{method = "d1"} chooses the derivative-free root-finding and \code{method = "d2"} chooses the
#' first-order root-finding.
#' @author Xiangjie Xue
solvegradientmultiple = function(x, mu0, pi0, points, tol = 1e-6, method = "auto"){
  # Here mu0 and pi0 is always fixed. since R has a copying sematics, we can modify
  # the x object so that it pre-computes flexden to save some time.
  # the corresponding methods in gradidentfunction should have capacity for testing whether
  # there is a precompute value
  # whether to precompute fullden remain a mystery
  if (!is.null(x$fn))
    x$flexden = x$fn(mu0, pi0)
  switch(method,
         "auto" = {
          # minimal sample to test implementation
          xx = makeobject(points[1], method = attr(x, "class"), beta = x$beta) # in the case the structure parameter is compulsory.
          rr = gradientfunction(xx, points[1], mu0 = 0, pi0 = 1, order = c(1, 1, 1))
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
  iter = 1;
  x = makeobject(x, pi0 = init, method = attr(x, "class"))
  r = computemixdist(x, mix = mix, maxiter = maxiter, tol = tol)
  d1 = r$ll + val
  d = estpi0d(x, r$mix$pt, r$mix$pr)
  d2 = d$d2; d3 = d$d3
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
      if ((init - neg) / (pos - neg) < 0.5){
        s = (pos + neg) / 2
        x = makeobject(x, pi0 = s, method = attr(x, "class"))
        r2 = computemixdist(x, mix = mix, maxiter = maxiter, tol = tol)

        if (r2$ll + val < 0){
          neg = s
        }else{
          pos = s
        }
      }
      neg = max(neg, init)
    }else{
      if ((pos - init) / (pos - neg) < 0.5){
        s = (pos + neg) / 2
        x = makeobject(x, pi0 = s, method = attr(x, "class"))
        r2 = computemixdist(x, mix = mix, maxiter = maxiter, tol = tol)

        if (r2$ll + val < 0){
          neg = s
        }else{
          pos = s
        }
      }
      pos = min(pos, init)
    }

    if (is.null(d3)){
      init = init - d1 * (1 - init) / d2;
    }else{
      init = init - 2 * d1 * d2 * (1 - init) / (2 * d2^2 - d1 * d3)
    }


    if (init < neg | init > pos){
      init = (pos + neg) / 2
    }

    x = makeobject(x, pi0 = init, method = attr(x, "class"))
    j0 = r$mix$pt == 0
    newmix = list(pt = r$mix$pt[!j0], pr = r$mix$pr[!j0] / sum(r$mix$pr[!j0]))
    r = computemixdist(x, mix = newmix, maxiter = maxiter, tol = tol)
    d1 = r$ll + val
    d = estpi0d(x, r$mix$pt, r$mix$pr)
    d2 = d$d2; d3 = d$d3
  }

  r
}
