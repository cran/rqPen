\name{scad_2_deriv}
\alias{scad_2_deriv}
\title{SCAD 2 Deriv}
\usage{
scad_2_deriv(x,lambda,a)
}
\arguments{
  \item{x}{Number to be evaluated.}
  \item{lambda}{Tuning parameter lambda.}
  \item{a}{Tuning paramter a. Default value of 3.7 as suggested in Fan and Li (2001).}
}
\value{
SCAD function is the difference of two convex functions. This is the derivative of the negative portion.
}
\description{
Derivative of the negative portion of the SCAD penalty function as described in Fan and Li (2001).
}
\examples{
scad_2_deriv(3,1)
scad_2_deriv(-3,1)
scad_2_deriv(.001,2)
}
\references{
[1] Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties, \emph{J. Am. Statist. Ass}, \bold{96}, 1348--1360.

[2] Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.
}
\author{Ben Sherwood}