\name{lilliefors.exp}
\alias{lilliefors.exp}
\alias{psi.plks.exp}
\alias{psi.qlks.exp}
\alias{psi.rlks.exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Approximate distribution of the Lilliefors's statistic for exponential
  hypotheses
}
\description{
  Distribution function, quantile function and random
  generation for the distribution of the Lilliefors's statistic for
  exponential hypotheses.
}
\usage{
psi.plks.exp(q, df, Nsim = 1000, lower.tail = TRUE)
psi.qlks.exp(p, df, Nsim = 1000)
psi.rlks.exp(n, df)
}
\arguments{
  \item{p}{  vector of probabilities. }
  \item{q}{ vector of quantiles. }
  \item{df}{ degrees of freedom of the statistic.  }
  \item{n}{ number of observations. }
  \item{Nsim}{ number of simulations. }
  \item{lower.tail}{ if TRUE (default) lower tail, otherwise upper tail. }
}
\details{
  The prefixes 'p', 'q' and 'r' follow the conventions of distribution
  functions in base \R. The prefix 'psi.' is a reminder that they are
  in package "psistat".

  The distribution function and the quantiles are obtained from a sample
  of size \code{Nsim}, ideally large. So, they should not be expected to
  give identical results in repeated calls with the same parameters.

  Instances of the statistic are obtained by generating random samples
  of length \code{df} from the exponential distribution, transforming
  each by the formula \eqn{z=x/{\bar x}}{z=x/xbar}, and calculating the
  Kolmogorov statistic for the transformed sample, \eqn{z}, with
  \eqn{F_0}{F0} set to the exponential distribution with rate 1.

  todo: pdf
}
\value{
  for \code{psi.plks.exp}, a numeric vector of the same length as the
  vector \code{q}.

  for \code{psi.qlks.exp}, a numeric vector of the same length as the
  vector \code{p}.

  for \code{psi.rlks.exp}, a numeric vector of length \code{n}
  representing a sample from the distribution of the Lilliefors KS
  statistic for the exponential distribution.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Georgi N. Boshnakov}
\seealso{
  \code{lillie.test} in package \code{nortest} for the Lilliefors test
  for normality.
}
\examples{
psi.plks.exp(c(0.1,0.2,0.3,0.4), 10, Nsim=1000)
psi.qlks.exp(c(0.25,0.5,0.75), 10, Nsim=10000)

xi <- seq(0,1,length=100)
yi <- psi.qlks.exp(xi,10,Nsim=10000)
plot(xi,yi)


# aproximate cdf of Lilliefors stat.
plot(ecdf(psi.rlks.exp(1000,10)))

# overlay the cdf of Kolmogorov for comparison.
lines(seq(0,1,length=100), psi.pks(seq(0,1,length=100),10),col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
