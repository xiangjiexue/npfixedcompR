# This file contains general functions used in npnorm family
# functions for normal mixture

dnpnorm = function(x, mu0 = 0, pi0 = 1, sd = 1, log = FALSE){
  # This version explicitly allow subprobability measures.
  # Hence no check on pi0
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = .rowSums(dnorm(x, mean = rep(mu0, rep(length(x), length(mu0))), sd = sd) *
                    rep(pi0, rep(length(x), length(pi0))), m = length(x), n = length(mu0))
  if (log) log(temp) else temp
}

pnpnorm = function(x, mu0 = 0, pi0 = 1, sd = 1, lower.tail = TRUE, log.p = FALSE){
  # This version explicitly allow subprobability measures.
  # Hence no check on pi0
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  temp = .rowSums(pnorm(x, mean = rep(mu0, rep(length(x), length(mu0))), sd = sd, lower.tail = lower.tail) *
                    rep(pi0, rep(length(x), length(pi0))), m = length(x), n = length(mu0))
  if (log.p) log(temp) else temp
}

pnpnorm1 = function(x, mu0 = 0, pi0 = 1, sd = 1, lower.tail = TRUE, log.p = FALSE){
  if (length(mu0) != length(pi0))
    stop("Length mismatch")
  j0 = pi0 == 0
  if (sum(!j0) > 0){
    temp = pnorm(x, mean = rep(mu0[!j0], rep(length(x), sum(!j0))), sd = sd, lower.tail = lower.tail, log.p = TRUE) +
      rep(log(pi0[!j0]), rep(length(x), sum(!j0)))
    dim(temp) = c(length(x), sum(!j0))
    maxcoef = apply(temp, 1, max)
    temp = log(rowSums(exp(temp - maxcoef))) + maxcoef
  }else{
    temp = rep(-Inf, length(x))
  }

  if (log.p) temp else exp(temp)
}

# functions of discrete normal mixture

logspace.add = function(lx, ly){
  j0 = lx == -Inf & ly == -Inf
  ans = pmax(lx, ly) + log1p(exp(-abs(lx - ly)))
  ans[j0] = -Inf
  ans
}

log1mexp = function(x){
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

logspace.sub = function(lx, ly){
  lx + log1mexp(lx - ly)
}

ddiscnorm = function(x, mean = 0, sd = 1, h = 1, log = FALSE){
  LLL = max(length(x), length(mean))
  xx = rep(x, length.out = LLL)
  meanx = rep(mean, length.out = LLL)
  temp = logspace.sub(pnorm(xx + h, meanx, sd, log.p = TRUE), pnorm(xx, meanx, sd, log.p = TRUE))
  if (log) temp else exp(temp)
}

dnpdiscnorm = function(x, mu0 = 0, pi0 = 0, sd = 1, h = 1, lower.tail = TRUE, log.p = FALSE){
  j0 = pi0 == 0
  if (sum(!j0) > 0){
    temp = logspace.sub(pnorm(x + h, rep(mu0[!j0], rep(length(x), sum(!j0))), sd = sd, lower.tail = lower.tail, log.p = TRUE),
                        pnorm(x, rep(mu0[!j0], rep(length(x), sum(!j0))), sd = sd, lower.tail = lower.tail, log.p = TRUE)) +
      rep(log(pi0[!j0]), rep(length(x), sum(!j0)))
    dim(temp) = c(length(x), sum(!j0))
    maxcoef = apply(temp, 1, max)
    temp = log(rowSums(exp(temp - maxcoef))) + maxcoef
  }else{
    temp = rep(-Inf, length(x))
  }

  if (log.p) temp else exp(temp)
}

pdiscnorm = function(x, mean =  0, sd = 1, h = 1, lower.tail = TRUE, log.p = FALSE){
  pnorm(x + h, mean = mean, sd = sd, lower.tail = lower.tail, log.p = log.p)
}

pnpdiscnorm = function(x, mu0 = 0, pi0 = 0, sd = 1, h = 1, lower.tail = TRUE, log.p = FALSE){
  pnpnorm1(x + h, mu0 = mu0, pi0 = pi0, sd = sd, lower.tail = lower.tail, log.p = log.p)
}

gridpointsnpnorm = function(x, grid=100) {
  rx = range(x$v)
  breaks = pmax(ceiling(diff(rx) / (5*x$beta)), 5)   # number of breaks
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
