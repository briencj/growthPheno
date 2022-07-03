\name{anom}
\alias{anom}
\title{Tests if any values in a vector are anomalous in being outside specified limits}
\description{Test whether any values in \code{x} are less than the value of 
             \code{lower}, if it is not \code{NULL}, or are greater than the 
             value of \code{upper}, if it is not \code{NULL}, or both.}
\usage{anom(x, lower=NULL, upper=NULL, na.rm = TRUE)}
\arguments{
 \item{x}{A \code{\link{vector}} containing the values to be tested.}
 \item{lower}{A \code{\link{numeric}} such that values in \code{x} below it are  
              considered to be anomalous.}
 \item{upper}{A \code{\link{numeric}} such that values in \code{x} above it are  
              considered to be anomalous.}
 \item{na.rm}{A \code{\link{logical}} indicating whether \code{NA} values should 
              be stripped before the testing proceeds.}
} 
\value{A \code{\link{logical}} indicating whether any values have been found to 
       be outside the limits specified by \code{lower} or \code{upper} or both.}
\author{Chris Brien}
\examples{
data(exampleData)
anom.val <- anom(longi.dat$sPSA.AGR, lower=2.5)
}
\keyword{data}
\keyword{manip}