\name{psi.rbern}
\alias{psi.rbern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate a random sample from Bernoulli distribution}
\description{Simulate a random sample from the Bernoulli distribution.}
\usage{
psi.rbern(n, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations.}
  \item{p}{probability of success.}
}
\details{
  A random variable has the Bernoulli(p) distribution if takes the
  value one with probability \eqn{p}  and the value zero with
  probability \eqn{1-p}.

  The Bernoulli(p) distribution is a particular case of the binomial
  distribution with \eqn{n=1}. 

}
\value{
  a numeric vector representing a random sample from Bernoulli(p)
}
\author{Georgi N. Boshnakov (but note that this function is a simple wrapper for
  \code{rbinom} from core \R)}
\examples{
psi.rbern(10, 0.5)   # 10 draws of a fair coin
psi.rbern(10, 0.6)   # 10 draws of a biased coin

# eqivalent to the above (but don't expect identical results!)
rbinom(10, size = 1, prob = 0.5)
rbinom(10, size = 1, prob = 0.6)
}
\keyword{distribution}
% __ONLY ONE__ keyword per line
