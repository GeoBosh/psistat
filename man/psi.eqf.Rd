\name{psi.eqf}
\alias{psi.eqf}
\title{ Empirical quantile function }
\description{
  Computes the empirical quantile function for a data vector and returns it
  as a step function object.
}
\usage{
psi.eqf(x, val.right = NULL)
}
\arguments{
  \item{x}{ a numeric data vector representing a random sample from a
    continuous distribution }
  \item{val.right}{ value of eqf outside [0,1], see details }
}
\details{
  The eqf for a sample of size n is a step function with value equal to
  the i-th order statistic in the interval ((i-1)/n,i/n] for i=1,...,n.

  Other definitions of sample quantiles, see \code{\link{quantile}},
  are not implemented by this function.
}
\value{
  The eqf as a \code{stepfun} object.
}
\author{Georgi N. Boshnakov}
\note{The result will be incorrect if \code{x} contains ties but this
  should not happen for continuous distributions.}
\seealso{
  \code{\link{quantile}},
  \code{\link{psi.plot.stepfun}},
  \code{\link{plot.stepfun}},
  \code{\link{stepfun}},
  \code{\link{ecdf}} }
\examples{
f5 <- psi.eqf(rnorm(20))
psi.plot.stepfun(f5,xlim=c(0,1),rigid.xlim=TRUE,verticals=FALSE)
psi.plot.stepfun(f5,xlim=c(0,1),rigid.xlim=TRUE,verticals=FALSE,pch=19)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
