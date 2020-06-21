dnpt = function(x, mu0 = 0, pi0 = 0, df, log = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = .rowSums(dt(x, ncp = rep(mu0, rep(length(x), length(mu0))), df = df) *
                    rep(pi0, rep(length(x), length(pi0))), m = length(x), n = length(mu0))
  if (log) log(temp) else temp
}

pnpt = function(x, mu0 = 0, pi0 = 0, df, lower.tail = TRUE, log.p = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = .rowSums(pt(x, ncp = rep(mu0, rep(length(x), length(mu0))), df = df, lower.tail = lower.tail) *
                    rep(pi0, rep(length(x), length(pi0))), m = length(x), n = length(mu0))
  if (log.p) log(temp) else temp
}

gridpointsnpt = function(x, grid=100) {
  rx = range(x$v)
  fac = ifelse(is.finite(x$beta), sqrt(x$beta / (x$beta - 2)), 1)
  breaks = pmax(ceiling(diff(rx) / (5 * fac)), 5)   # number of breaks
  if (is.null(x$w)) {w = rep(1, length(x$v))} else {w = x$w}
  r = whist(x$v, w, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
  i = r$density != 0
  i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
  m = sum(i)
  k = pmax(ceiling(grid / m), 10)           # at least 10 in each interval
  d = r$breaks[2] - r$breaks[1]
  s = r$breaks[-length(r$breaks)][i]
  sort(c(rx[1], rep(s, rep(k,m)) + d * (1:k-0.5)/k, rx[2]), decreasing = FALSE)
}

#' @rdname makeobject
#' @export
makeobject.nptll = function(v, mu0, pi0, beta){
  if (class(v) == "nptll"){
    # update information
    x = v
    if (!missing(mu0)) x$mu0 = mu0
    if (!missing(pi0)) x$pi0 = pi0
    if (!missing(beta)) x$beta = Inf
    x$precompute = dnpt(x$v, mu0 = x$mu0, pi0 = x$pi0, df = x$beta)
  }

  if (is.numeric(v)){
    if (missing(mu0)) mu0 = 0
    if (missing(pi0)) pi0 = 0
    if (missing(beta)) beta = Inf
    x = list(v = v, mu0 = mu0, pi0 = pi0, beta = beta,
             precompute = dnpt(v, mu0 = mu0, pi0 = pi0, df = beta))
    attr(x, "class") = "nptll"
  }

  x
}

lossfunction.nptll = function(x, mu0, pi0){
  -sum(log(dnpt(x$v, mu0 = mu0, pi0 = pi0, df = x$beta) + x$precompute))
}

gradientfunction.nptll = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  flexden = dnpt(x$v, mu0 = mu0, pi0 = pi0, df = x$beta)
  temp = dt(x$v, ncp = rep(mu, rep(length(x$v), length(mu))), df = x$beta) * sum(pi0)
  fullden = flexden + x$precompute
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")

  # only d0 is implemented
  if (order[1] == 1){
    ans$d0 = .colSums((flexden - temp) / fullden, m = length(x$v), n = length(mu))
  }


  ans
}

#' @export
computemixdist.nptll = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    fac = ifelse(is.finite(x$beta), sqrt(x$beta / (x$beta - 2)), 1)
    breaks = pmax(ceiling(diff(rx) / (5 * fac)), 10)   # number of breaks
    r = whist(x$v, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    r$density = pmax(0, r$density - pnpt(r$breaks[-1], mu0 = x$mu0, pi0 = x$pi0, df = x$beta) +
                       pnpt(r$breaks[-length(r$breaks)], mu0 = x$mu0, pi0 = x$pi0, df = x$beta))
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0] / sum(r$density)
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = gridpointsnpt(x)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol))
    pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
    sp = dt(x$v, ncp = rep(mu0new, rep(length(x$v), length(mu0new))), df = x$beta)
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
             family = "npt",
             max.gradient = min(gradientfunction(x, mu0, mu0, pi0, order = c(1, 0, 0))$d0),
             mix = list(pt = r$pt, pr = r$pr),
             ll = nloss,
             dd0 = gradientfunction(x, 0, r$pt, r$pr, order = c(1, 0, 0))$d0,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}

#' @rdname estpi0
#' @export
estpi0.nptll = function(x, val = 0.5 * log(length(x$v)), mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  x = makeobject(x, pi0 = 1 - tol / 2, method = attr(x, "class"))
  r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
  x = makeobject(x, pi0 = 0, method = attr(x, "class"))
  r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)

  if (r1$ll - r0$ll < val){
    r = list(iter = 0,
             family = "npt",
             max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             mix = list(pt = 0, pr = 1),
             ll = lossfunction(x, mu0 = 0, pi0 = 1),
             dd0 = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
             convergence = 0)
  }else{
    r = solveestpi0(x = x, init = dnpt(0, mu0 = r0$mix$pt, pi0 = r0$mix$pr, df = x$beta) / dt(0, df = x$beta),
                    val = -r0$ll - val, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
  }

  r
}
