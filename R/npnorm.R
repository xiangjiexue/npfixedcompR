# This file contains general functions used in npnorm family

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

# taken from nspmix with modification
gridpointsnpnorm = function(x, grid=100) {
  rx = range(x$v)
  breaks = pmax(ceiling(diff(rx) / (5*x$sd)), 5)   # number of breaks
  r = hist(x$v, breaks = breaks, probability = TRUE, plot = FALSE, warn.unused = FALSE)
  i = r$density != 0
  i = i | c(i[-1],FALSE) | c(FALSE,i[-length(i)])  # include neighbours
  m = sum(i)
  k = pmax(ceiling(grid / m), 10)           # at least 10 in each interval
  d = r$breaks[2] - r$breaks[1]
  s = r$breaks[-length(r$breaks)][i]
  sort(c(rx[1], rep(s, rep(k,m)) + d * (1:k-0.5)/k, rx[2]), decreasing = FALSE)
}
