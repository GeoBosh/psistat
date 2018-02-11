\name{psi.Dn}
\alias{psi.Dn}
\title{ Compute Kolmogorov statistics }
\description{
  Compute the Kolmogorov statistic and its one-sided variants.
}
\usage{
psi.Dn(x, cdf = pnorm, type = c("D!=", "D^+", "D^-"),
       ties.jitter = NULL, ...)
}
\arguments{
  \item{x}{ numeric; the data vector }
  \item{cdf}{ a character string naming a cumulative distribution
    function or an actual cumulative distribution function such as
    'pnorm'.}

  \item{type}{ specifies the required variant of the statistic and must
    be one of  '"\code{D!=}"' (default) for \eqn{D_n}{Dn}, '"\code{D^+}"' for
    \eqn{D^+}{D+}, or '"\code{D^-}"' for \eqn{D^-}{D-}.  }

  \item{ties.jitter}{ processing of  ties: if NULL (default) no check
    for ties is done, if TRUE a default jitter is applied to tied values,
    a numeric value, say \code{a}, results  in adding a U(-a,a) r.v. to
    the tied observations.
  }

  \item{\dots}{ parameters of the distribution function specified by
  \code{cdf}.}
}
\details{
  The  same result may be obtained by extracting the \code{statistic}
  component from the result of \code{ks.test}, see the examples.
  \code{psi.Dn} may be convenient for use as part of other calculations
  since it only computes \eqn{D}.

  When there are ties in the data and \code{ties.jitter} is not NULL the
  returned value is random and will vary even if the function is
  called with the same data.

}
\value{
  the value of the requested statistic, a single number.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Georgi N. Boshnakov}
%\note{ ~~further notes~~}
\seealso{ \code{\link{ks.test}}  }
\examples{
a1 <- rnorm(10)                # some data
psi.Dn(a1,cdf=pnorm)           # Compute Dn
ks.test(a1,pnorm)$statistic    # should  be the same

psi.Dn(a1,cdf=pnorm,type="D^+")
ks.test(a1,pnorm,alternative="greater")$statistic

psi.Dn(a1,cdf=pnorm,type="D^-")
ks.test(a1,pnorm,alternative="less")$statistic

psi.Dn(a1,cdf=pnorm,mean=2,sd=2.3)  # pnorm with parameters.


a2 <- rexp(10,rate=0.1)           # exponential dist
psi.Dn(a1,cdf=pexp)
psi.Dn(a1,cdf=pexp,rate=0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
