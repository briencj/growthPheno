\name{validSmoothsFrame}
\alias{validSmoothsFrame}
\title{Checks that an object is a valid \code{\link{smooths.frame}}.}
\description{Checks that an \code{object} is a \code{\link{smooths.frame}} 
             of S3-class \code{data.frame} that contains the columns 
             \code{Type}, \code{TunePar}, \code{TuneVal}, \code{Tuning}, 
             \code{Method}, as well as the columns specified by the 
             atttributes of the \code{object}, namely \code{individuals} 
             and \code{times}.}

\usage{
validSmoothsFrame(object)}
\arguments{
 \item{object}{a \code{\link{smooths.frame}}.}
}
\value{\code{TRUE} or a \code{character} describing why the \code{object} 
       is not a valid \code{\link{smooths.frame}}.}


\author{Chris Brien}
\seealso{\code{\link{is.smooths.frame}}, \code{\link{as.smooths.frame}}}
\examples{
dat <- read.table(header = TRUE, text = "
Type TunePar TuneVal Tuning Method       ID  DAP   PSA      sPSA
NCSS      df       4   df-4 direct 045451-C   28 57.446 51.18456
NCSS      df       4   df-4 direct 045451-C   30 89.306 87.67343
NCSS      df       7   df-7 direct 045451-C   28 57.446 57.01589
NCSS      df       7   df-7 direct 045451-C   30 89.306 87.01316
")
dat[1:7] <- lapply(dat[1:6], factor)
dat <- as.smooths.frame(dat, individuals = "ID", times = "DAP")
is.smooths.frame(dat)
validSmoothsFrame(dat)
}
\keyword{asreml}
\keyword{htest}