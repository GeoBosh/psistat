\name{psi.jitter}
\alias{psi.jitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Jitter ties }
\description{
  Apply jitter to the tied observations only.
}
\usage{
psi.jitter(x, amount)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ numeric vector of observations. }
  \item{amount}{ amount of jitter, see details. }
}
\details{
  Add random numbers U(-amount,amount) to the tied  observations.
}
\value{
  a numeric vector.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Georgi N. Boshnakov}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
