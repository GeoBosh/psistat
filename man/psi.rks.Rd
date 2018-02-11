\name{psi.rks}
\alias{psi.rks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random samples from Kolmogorov distribution }
\description{
  Generate random samples from Kolmogorov distribution.
}
\usage{
psi.rks(n, df, randgen = runif, cdf = punif, ...)
}
\arguments{
  \item{n}{ number of observations. }
  \item{df}{ degrees of freedom of the Kolmogorov distribution  }
  \item{randgen}{a function or a name of one that generates random
    numbers from the  hypothesised distribution. }
  \item{cdf}{ the cdf of the hypothesised distribution. }
  \item{\dots}{ parameters to be passed down to \code{randgen} and \code{cdf} }
}
\details{
  Instances of the \eqn{D_{df}} statistic are obtained by generating
  random samples of length \code{df} from  the hypothesised distribution
  and computing \eqn{D_{df}} on them.

  In the case of a simple null hypothesis from a continuous
  distribution the distribution of the statistic  does not depend on
  the underlying distribution. So, only the arguments \code{n} and
  \code{df} are really necessary for  the Kolmogorov-Smirnov test.
  Imaginative settings for this argument  may be used to generate samples
  in the case of  composite  hypotheses.
}
\value{
  a  numeric vector representing a sample from the distribution of the
  Kolmogorov statistic for the specified distribution.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Georgi N. Boshnakov}
%\note{ ~~further notes~~}
\seealso{ \code{\link{psi.pks}}%
  %, \code{\link{psi.rlks}}%
  %, etc.
}
\examples{
x1 <- psi.rks(1000,10)
hist(x1)

plot(ecdf(x1))                   # simulation estimate of the dist of D_10

f1 <- function(x) psi.pks(x,10)
curve(f1,0,1,col="blue",add=TRUE)  # overlay the exact cdf
}
\keyword{ distribution }
