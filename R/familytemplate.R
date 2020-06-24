# This file contains the template for implementing a new family. Make changes accordingly.
# If a new density is used rather than npnorm and npt, you need to define it.
# This file uses the term new-family as a new method, change when implementing (can use find all & replace)

# makeobject.newfamily = function(v, mu0, pi0, beta){
#   if (class(v) == "newfamily"){
#     x = v
#     if (!missing(mu0)) x$mu0 = mu0
#     if (!missing(pi0)) x$pi0 = pi0
#     if (!missing(beta)) x$beta = beta
#     # updating information in v should be included below, using x$...
#   }
#
#   if (is.numeric(v)){
#     if (missing(mu0)) mu0 = 0
#     if (missing(pi0)) pi0 = 0
#     if (missing(beta)) beta = 1
#     # construction. Include any precomputed value as ...
#     x = list(v = v, mu0 = mu0, pi0 = pi0, beta = beta, ...)
#     attr(x, "class") = "newfamily"
#   }
#
#   x
# }

# lossfunction.newfamily = function(x, mu0, pi0){
#   # The loss function for minimisation
# }

# gradientfunction.newfamily = function(x, mu, mu0, pi0, order = c(1, 0, 0)){
#   # This is general for performance reason
#   if (!is.null(x$flexden)){
#     flexden = x$flexden
#   }else{
#     # flexden = any values you do not want to repeatedly compute.
#   }
#   ans = vector("list", 3)
#   names(ans) = c("d0", "d1", "d2")
#   if (order[1] == 1){
#     # ans$d0 = expression for the gradient function (vectorised)
#   }
#   if (order[2] == 1){
#     # ans$d1 = expression for the first derivative of the gradient function (vectorised) (if needed)
#   }
#   if (order[3] == 1){
#     # ans$d2 = expression for the second derivative of the gradient function (vectorised) (if needed)
#   }
#
#   ans
# }

# computemixdist.newfamily = function(x, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
#   if (is.null(mix)){
#     # simple estimate of initial distribution if mix is not provided
#   }else{
#     mu0 = mix$pt; pi0 = mix$pr
#   }
#
#   # x$fn = the functions that computes any values you do not want to repeatedly compute, linking to gradient function
#
#   pi0 = pi0 * (1 - sum(x$pi0))
#
#   iter = 0; convergence = 0
#   closs = lossfunction(x, mu0, pi0)
#   # points = the grid for searching new support points
#
#   repeat{
#     mu0new = c(mu0, solvegradientmultiple(x, mu0, pi0, points, tol, method = "d1")) # method can be chosen if you know which is most efficient, or use "auto"
#     pi0new = c(pi0, rep(0, length.out = length(mu0new) - length(mu0)))
#
#     # chunk for estimating proportion vectors based on the new support set.
#
#     # use checklossfunction for line search if the second order expansion is not exact.
#
#     r = collapsemix(x, r$pt, r$pr, tol) # collapsing
#     mu0 = r$pt; pi0 = r$pr
#     iter = iter + 1
#     nloss = lossfunction(x, mu0, pi0)
#
#     if (verbose){
#       cat("Iteration: ", iter, "\n")
#       cat(paste0("Support Point ", round(r$pt, -ceiling(log10(tol))), " with probability ", round(r$pr, -ceiling(log10(tol))), "\n"))
#       cat("Current log-likelihood ", as.character(round(-nloss, -ceiling(log10(tol)))), "\n")
#     }
#
#     if (closs - nloss < tol){
#       convergence = 0
#       break
#     }
#
#     if (iter > maxiter){
#       convergence = 1
#       break
#     }
#     closs = nloss
#   }
#
#   mu0new = c(mu0, x$mu0); pi0new = c(pi0, x$pi0)
#   index = order(mu0new, decreasing = FALSE)
#   r = unique.disc(mu0new[index], pi0new[index])
#
#   ans = list(iter = iter,
#              family = "npnorm", # specifying family. useful in implementing rejectregion and posterior mean
#              max.gradient = min(gradientfunction(x, mu0, mu0, pi0, order = c(1, 0, 0))$d0),
#              mix = list(pt = r$pt, pr = r$pr),
#              ll = nloss,
#              beta = x$beta,
#              dd0 = gradientfunction(x, 0, mu0, pi0, order = c(1, 0, 0))$d0,
#              convergence = convergence)
#
#   attr(ans, "class") = "nspmix"
#   ans
# }

# estimating the null proportion (optional)
# estpi0.newfamily = function(x, val, mix = NULL, tol = 1e-6, maxiter = 100, verbose = FALSE){
#   x = makeobject(x, pi0 = 1 - tol / 2, method = attr(x, "class"))
#   r1 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
#   x = makeobject(x, pi0 = 0, method = attr(x, "class"))
#   r0 = computemixdist(x, mix = mix, tol = tol, maxiter = maxiter)
#
#   if (condition){
#     r = list(iter = 0,
#              family = "npnorm",# specifying family. useful in implementing rejectregion and posterior mean
#              max.gradient = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
#              mix = list(pt = 0, pr = 1),
#              beta = x$beta,
#              ll = lossfunction(x, mu0 = 0, pi0 = 1),
#              dd0 = gradientfunction(x, 0, 0, 1, order = c(1, 0, 0))$d0,
#              convergence = 0)
#   }else{
#     r = solveestpi0(x = x, init = initialvalue,
#                     val = targetvalueforfindingzero, mix = r0$mix, tol = tol, maxiter = maxiter, verbose = verbose)
#   }
#
#   r
# }
