\name{psi.pks}
\alias{psi.pks}
\title{ The Kolmogorov Distribution }
\description{
  Compute the cdf of the Kolmogorov statistic for the two-sided test.
}
\usage{
psi.pks(q, n, lower.tail = TRUE, log.p = FALSE, exact = NULL, warn = FALSE)
}
\arguments{
  \item{q}{  vector of quantiles. }
  \item{n}{ sample size. }
  \item{lower.tail}{ logical; if TRUE (default) probabilities are P[X <= x],
          otherwise, P[X > x]. }
  \item{log.p}{ logical; if TRUE, probabilities p are given as log(p). }
  \item{exact}{if TRUE use exact method, if FALSE use asymptotic
    aproximation, if NULL (default) use exact method for n<=100 and
    asymptotic for n>100.}
  \item{warn}{ logical; if TRUE (default) and \code{exact=TRUE} and \code{n>100}
    print a message.}
}
\details{
  Avoid setting exact=TRUE for large values of n
  since the exact method can become very slow for such n.

  This function is a wrapper for two undocumented low level \code{R}
  functions which are used by \code{ks.test}.
}
\value{
  numerical vector of same length as \code{q} containing the values of
  the cdf of Kolmogorov's distribution.
}
\references{ See \code{ks.test} and, in the R sources, \code{ks.c} }
\author{Georgi N. Boshnakov; the essential computation is done by C
  functions by R Core Development Team in package stats (copied here
  since unexported).}
% \note{}
\seealso{ \code{\link{ks.test}} }
\examples{
xi <- seq(0,1,length=100)            # some x values

plot(xi,psi.pks(xi,4))               # cdf of D_4
lines(xi,psi.pks(xi,4))
lines(xi,psi.pks(xi,50),col="blue")  # overlay the cdf of  D_{50}
lines(xi,psi.pks(xi,100),col="red")  # overlay the cdf of  D_{100}
abline(h=0.95, col="brown")

# q : P(Dn < q) = 0.95 for selected values of n
psi.pks(0.6239385,4)
psi.pks(0.2940753,20)
psi.pks(0.1340279,100)
psi.pks(0.04294685,1000)   # asymptotic approximation


n <- 8
plot(xi,psi.pks(xi,n,exact=FALSE))    # evaluate and plot asymptotic approx.
lines(xi,psi.pks(xi,n))               # overlay the exact cdf.

f1 <- function(x) psi.pks(x,10)      # f1(x) computes cdf of D_10
f2 <- function(x) psi.pks(x,20)      # f2(x) computes cdf of D_20
f3 <- function(x) psi.pks(x,100)     # f3(x) computes cdf of D_100

curve(f1,0,1)
curve(f2,0,1,add=TRUE,col="blue")
curve(f3,0,1,add=TRUE,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
