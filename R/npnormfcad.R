#' @rdname makeobject
makeobject.npnormad = function(v, mu0 = 0, pi0 = 0, beta = 1){
  if (class(v) == "npnormad"){
    # update information
    v$mu0 = mu0; v$pi0 = pi0; v$sd = beta
    return(v)
  }

  if (is.numeric(v)){
    x = list(v = v, mu0 = mu0, pi0 = pi0, sd = beta,
             a1 = seq(from = 1 / length(v), to = (2 * length(v) - 1) / length(v), length = length(v)),
             a2 = seq(from = (2 * length(v) - 1) / length(v), to = 1 / length(v), length = length(v)))
    class(x) = "npnormad"
    return(x)
  }
}

lossfunction.npnormad = function(x, mu0, pi0){
  temp = pnpnorm(x$v, mu0, pi0, x$sd)
  -sum(log(temp) * x$a1 + log1p(-temp) * x$a2)
}
