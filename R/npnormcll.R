dnormc = function(x, mean = 0, n, log = FALSE){
  if (any(abs(x) >= 1) | any(abs(mean) >= 1))
    stop("Error in specifying data or mean")
  LLL = max(length(x), length(mean))
  xx = rep(x, length.out = LLL)
  meanx = rep(mean, length.out = LLL)
  dnorm(xx, mean = meanx, sd = (1 - meanx^2) / sqrt(n), log = log)
}

pnormc = function(x, mean = 0, n, lower.tail = TRUE, log.p = FALSE){
  if (any(abs(x) >= 1) | any(abs(mean) >= 1))
    stop("Error in specifying data or mean")
  LLL = max(length(x), length(mean))
  xx = rep(x, length.out = LLL)
  meanx = rep(mean, length.out = LLL)
  pnorm(xx, mean = meanx, sd = (1 - meanx^2) / sqrt(n), lower.tail = lower.tail, log.p = log.p)
}

dnpnormc = function(x, mu0 = 0, pi0 = 0, n, log = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = dnormc(x, mean = rep(mu0, rep(length(x), length(mu0))), n = n)
  dim(temp) = c(length(x), length(mu0))
  temp = drop(temp %*% pi0)
  if (log) log(temp) else temp
}

pnpnormc = function(x, mu0 = 0, pi0 = 0, n, lower.tail = TRUE, log.p = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = pnormc(x, mean = rep(mu0, rep(length(x), length(mu0))), n = n, lower.tail = lower.tail)
  dim(temp) = c(length(x), length(mu0))
  temp = drop(temp %*% pi0)
  if (log.p) log(temp) else temp
}

# here beta is n
#' @rdname makeobject
#' @export
makeobject.npnormcll = function(v, mu0, pi0, beta){
  if (class(v) == "npnormcll"){
    x = v
    if (!missing(mu0)) x$mu0 = mu0
    if (!missing(pi0)) x$pi0 = pi0
    if (!missing(beta)) x$beta = beta
    # updating information in v should be included below, using x$...
    x$precompute = dnpnormc(x = x$v, mu0 = x$mu0, pi0 = x$pi0, n = x$beta)
  }

  if (is.numeric(v)){
    if (missing(mu0)) mu0 = 0
    if (missing(pi0)) pi0 = 0
    # if (missing(beta)) beta = 1
    # construction. Include any precomputed value as ...
    x = list(v = v, mu0 = mu0, pi0 = pi0, beta = beta,
             precompute = dnpnormc(x = v, mu0 = mu0, pi0 = pi0, n = beta))
    attr(x, "class") = "npnormcll"
  }

  x
}

lossfunction.npnormcll = function(x, mu0, pi0){
  -sum(log(dnpnormc(x$v, mu0 = mu0, pi0 = pi0, n = x$beta) + x$precompute))
}

gradientfunction.npnormcll = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
  # This is general for performance reason
  if (!is.null(x$flexden)){
    flexden = x$flexden
  }else{
    flexden = dnpnormc(x = x$v, mu0 = mu0, pi0 = pi0, n = x$beta)
  }
  murep = rep(mu, rep(length(x$v), length(mu)))
  temp = dnormc(x$v, mean = murep, n = x$beta) * sum(pi0)
  fullden = flexden + x$precompute
  ans = vector("list", 3)
  names(ans) = c("d0", "d1", "d2")
  if (order[1] == 1){
    ans$d0 = .colSums((flexden - temp) / fullden, m = length(x$v), n = length(mu))
  }
  if (order[2] == 1){
    temp2 = ((x$beta + 4) * murep^3 - 3 * x$beta * murep^2 * x$v + murep * (2 * x$beta * x$v^2 + x$beta - 2) -
               x$beta * x$v - 2 * murep^5) / (1 - murep^2)^3
    ans$d1 = .colSums(temp2 * temp / fullden, m = length(x$v), n = length(mu))
  }

  ans
}

#' @export
computemixdist.npnormcll = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
  if (is.null(mix)){
    rx = range(x$v)
    breaks = pmax(ceiling(diff(rx) / (5 / sqrt(x$beta))), 10)   # number of breaks
    r = whist(x$v, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
    r$density = pmax(0, r$density  / sum(r$density) - pnpnormc(r$breaks[-1], mu0 = x$mu0, pi0 = x$pi0, n = x$beta) +
                       pnpnormc(r$breaks[-length(r$breaks)], mu0 = x$mu0, pi0 = x$pi0, n = x$beta))
    mu0 = r$mids[r$density != 0]
    pi0 = r$density[r$density != 0] / sum(r$density)
  }else{
    mu0 = mix$pt; pi0 = mix$pr
  }

  x$fn = function(mu0, pi0) dnpnormc(x$v, mu0, pi0, x$beta)
  pi0 = pi0 * (1 - sum(x$pi0))

  iter = 0; convergence = 0
  closs = lossfunction(x, mu0, pi0)
  points = gridpointsnpnorm(x)
  points = points[points > -1 & points < 1] # because the support space is (-1, 1)

  repeat{
    mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol))
    pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
    sp = dnormc(x$v, mean = rep(mu0new, rep(length(x$v), length(mu0new))), n = x$beta)
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
             family = "npnormc",
             max.gradient = min(gradientfunction(x, mu0, mu0, pi0, order = c(1, 0, 0))$d0),
             mix = list(pt = r$pt, pr = r$pr),
             ll = nloss,
             beta = x$beta,
             convergence = convergence)

  attr(ans, "class") = "nspmix"
  ans
}
