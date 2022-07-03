\name{is.smooths.frame}
\alias{is.smooths.frame}
\title{Tests whether an object is of class smooths.frame}
\description{A single-line \code{function} that tests whether an object is of \code{class} \code{\link{smooths.frame}}.}
\usage{is.smooths.frame(object)}
\arguments{
 \item{object}{An \code{object} to be tested.}
}
\value{A \code{logical}.}
\author{Chris Brien}
\seealso{\code{\link{validSmoothsFrame}}, \code{\link{as.smooths.frame}}}
\examples{
dat <- read.table(header = TRUE, text = "
Type TunePar TuneVal Tuning Method       ID  DAP   PSA      sPSA
NCSS      df       4   df-4 direct 045451-C   28 57.446 51.18456
NCSS      df       4   df-4 direct 045451-C   30 89.306 87.67343
NCSS      df       7   df-7 direct 045451-C   28 57.446 57.01589
NCSS      df       7   df-7 direct 045451-C   30 89.306 87.01316
")
dat[1:7] <- lapply(dat[1:7], factor)
dat <- as.smooths.frame(dat, individuals = "ID", times = "DAP")
is.smooths.frame(dat)
validSmoothsFrame(dat)
}
\keyword{manip}