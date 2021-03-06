% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npfixedcompR.R
\name{estpi0}
\alias{estpi0}
\title{Computing non-parametric mixing distribution with estimated proportion at 0}
\usage{
estpi0(x, ...)
}
\arguments{
\item{x}{a object from implemented family}

\item{...}{parameters above passed to the specific method.}
}
\description{
computing non-parametric mixing distribution with estimated proportion at 0
}
\details{
This is a function for computing non-parametric mixing
distribution with estimated proportion at 0. Different families will
have different threshold values.

The parameters are listed as follows:
\itemize{
\item val: Threshold value
\item mix: The initial proper mixing distribution.
\item tol: tolerance to stop the code.
\item maxiter: maximum iterations allowed.
\item verbose: logical; Whether to print the intermediate results.
}

It is not shown in the parameter section since various method have different
default thresold values and this function essentially calls the class method
in the object.

The full list of implemented families is in \code{\link{makeobject}}.
}
\examples{
data = rnorm(500, c(0, 2))
pi0 = 0.5
x = makeobject(data, method = "npnormll")
estpi0(x)
x = makeobject(data, method = "npnormllw")
estpi0(x)
x = makeobject(data, method = "npnormcvm")
estpi0(x)
x = makeobject(data, method = "npnormcvm")
estpi0(x)
x = makeobject(data, method = "npnormad")
estpi0(x)
x = makeobject(data, method = "nptll")
estpi0(x)
}
