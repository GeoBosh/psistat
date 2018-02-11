\name{psi.lks.exp.test}
\alias{psi.lks.exp.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lilliefors test for exponentiality }
\description{
  Performs Lilliefors test for exponentiality.
}
\usage{
psi.lks.exp.test(x, Nsim = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of observations. }
  \item{Nsim}{ number of simulations to perform, see details. }
  \item{\dots}{ other parameters. }
}
\details{
  \code{psi.lks.exp.test} performs a Lilliefors test for
  exponentiality, i.e. the null hypothesis is that \code{x} is a random
  sample from the exponential distribution with some (unspecified)
  degrees of freedom.

  The p-value is calculated through simulation controlled by the
  argument \code{Nsim} (larger \code{Nsim} should give a more reliable
  value).

  todo: remove the '...' argument, it should not be used.
}
\value{
   A list with class '"htest"' (modelled after the result of
  \code{\link{ks.test}}) containing the following components:

   \item{statistic}{ the value of the test statistic. }

   \item{p.value}{ the p-value of the test. }

   \item{alternative}{ a character string describing the alternative
     hypothesis. }

   \item{method}{ a character string indicating what type of test was
     performed. }

   \item{data.name}{ a character string giving the name(s) of the
     data. }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Georgi N. Boshnakov}
\note{

  Since the p-value is calculated through simulation, calling
  \code{psi.lks.exp.test} twice with the same data will not give
  identical results. The computed p-values should be only slightly
  different (especially with large \code{Nsim}) and will hardly ever
  lead to contradicting decisions.

  If you want the same result each time set the random seed before
  calling the function to a particular value, say \code{set.seed(1234)}.
}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
psi.lks.exp.test(rexp(10),Nsim=10000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
