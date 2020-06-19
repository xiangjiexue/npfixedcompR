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
    sapply(index, function(ddd){
      solvegradientsingle(x = x, mu0 = mu0, pi0 = pi0, lower = points[ddd], upper = points[ddd + 1], tol)
    })
  }else{NULL}
}
