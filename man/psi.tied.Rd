\name{psi.tied}
\alias{psi.tied}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identify tied observations }
\description{
  Determines which elements of a vector are tied.
}
\usage{
psi.tied(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector }
}
\details{
  Unlike \code{duplicated}, this function returns all tied elements.
}
\value{
  a logical vector with the value TRUE at the  positions of the tied elements.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Georgi N. Boshnakov}
%\note{ ~~further notes~~ }
\seealso{ \code{\link{duplicated}} }
\examples{
v1 <- sample(1:6,12,replace=TRUE)
psi.tied(v1)
duplicated(v1)
d1 <- data.frame(v=v1, tied=psi.tied(v1), dup=duplicated(v1))
d1$flag <- ifelse(d1$tied==d1$dup, "", "different")
d1    # show the difference between `psi.tied' and `duplicated'
}
\keyword{ logic }
