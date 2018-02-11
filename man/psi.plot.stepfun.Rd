\name{psi.plot.stepfun}
\alias{psi.plot.stepfun}
\title{ Plot piecewise constant functions }
\description{
  Plot a stepfun object.  The function is a modification of
  plot.stepfun from package stats but has an additional argument to
  force the domain of x-values rigidly.
}
\usage{
psi.plot.stepfun(x, xval, xlim, ylim = range(c(y, Fn.kn)),
    xlab = "x", ylab = "f(x)", main = NULL,
    add = FALSE, verticals = TRUE, do.points = TRUE,
    pch = par("pch"), col.points = par("col"), cex.points = par("cex"),
    col.hor = par("col"), col.vert = par("col"),
    lty = par("lty"), lwd = par("lwd"),
    rigid.xlim = FALSE, ...)
}
\arguments{
  \item{x}{ an R object inheriting from '"stepfun"', usually created by
     \code{stepfun} }

  \item{xlim}{ numeric(2), range of 'x' values to use, has sensible defaults. }


  \item{add}{  logical; if 'TRUE' only \emph{add} to an existing plot. }
  \item{verticals}{  logical;  if 'TRUE', draw vertical lines at steps. }

  \item{rigid.xlim}{ If \code{TRUE} respect \code{xlim} rigidly, see details. }

  \item{xval}{ see help page of \code{plot.stepfun}. }
  \item{ylim}{ see help page of \code{plot.stepfun}. }
  \item{xlab}{ see help page of \code{plot.stepfun}. }
  \item{ylab}{ see help page of \code{plot.stepfun}. }
  \item{main}{ see help page of \code{plot.stepfun}. }
  \item{do.points}{ see help page of \code{plot.stepfun}. }
  \item{pch}{ see help page of \code{plot.stepfun}. }
  \item{col.points}{ see help page of \code{plot.stepfun}. }
  \item{cex.points}{ see help page of \code{plot.stepfun}. }
  \item{col.hor}{ see help page of \code{plot.stepfun}. }
  \item{col.vert}{ see help page of \code{plot.stepfun}. }
  \item{lty}{ see help page of \code{plot.stepfun}. }
  \item{lwd}{ see help page of \code{plot.stepfun}. }
  \item{\dots}{ see help page of \code{plot.stepfun}. }
}
\details{
  The default method for plotting \code{stepfun} objects extends
  slightly the domain requested by \code{xlim}. This is not always
  desirable, especially if the function is not defined outside the
  specified limits. This function has all the arguments of
  \code{plot.stepfun} and does the same job with the additional option
  to force the use of \code{xlim} as given by setting the argument
  \code{rigid.xlim} to \code{TRUE}.
}
\value{
  A list with two components
  \item{t}{abscissa (x) values, including the two outermost ones.}
  \item{y}{y values `in between' the `t[]'.}
}
\references{ R package "stats" for the code of the original
  \code{plot.stepfun}.}  
\author{Georgi N. Boshnakov (to be blamed for bugs; the credits should
  go to the \R core team)}
\note{%
  This function is a modification of \code{plot.stepfun} from the
  \code{stats} package as supplied with R~2.8.1. Some of the text in
  this help page has been taken from the help page of
  \code{plot.stepfun}.
}
\seealso{\code{\link{plot.stepfun}}, \code{\link{stepfun}}}
\examples{
# define  empirical quantile functon as a step function.
eqf <- function(x) stepfun((1:(length(x)))/length(x),c(x,NA),right=TRUE)

# create eqf for a random sample.
x <- sort(rnorm(10))
f1 <- eqf(x)

# plot f1
psi.plot.stepfun(f1,xlim=c(0,1),rigid.xlim=TRUE)
psi.plot.stepfun(f1,xlim=c(0,1),rigid.xlim=TRUE,verticals=FALSE)
psi.plot.stepfun(f1,xlim=c(0,1),rigid.xlim=TRUE,verticals=FALSE,
                 main="An emprirical qf")

psi.plot.stepfun(f1,xlim=c(0,1),rigid.xlim=TRUE,pch=19,verticals=FALSE)

# plot(f1) would give an error because of the NA, so modify.
eqf2 <- function(x) stepfun((1:(length(x)))/length(x),c(x,0),right=TRUE)
f2 <- eqf2(x)

# the default method for stepfun plots outside the domain of eqf.
plot(f2,xlim=c(0,1),verticals=FALSE)

# eqf's with overlaid qf's
psi.plot.stepfun(f1,xlim=c(0,1),rigid.xlim=TRUE,pch=19,verticals=FALSE)
curve(qnorm,add=TRUE, col="red")

psi.plot.stepfun(eqf(sort(rnorm(100))),xlim=c(0,1),rigid.xlim=TRUE,
                 pch=19,verticals=FALSE,do.points=FALSE)
curve(qnorm,add=TRUE, col="red")
}
\keyword{ hplot }
