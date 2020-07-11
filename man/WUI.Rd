\name{WUI}
\alias{WUI}
\title{Calculates the Water Use Index (WUI)}
\description{Calculates the Water Use Index, returning \code{NA} if the water use is zero.}
\usage{WUI(response, water)}
\arguments{
 \item{response}{A \code{\link{numeric}} giving the value of the response achieved.}
 \item{water}{A \code{\link{numeric}} giving the amount of water used.}
} 
\value{A \code{\link{numeric}} containing the water divided by the response, unless water is 
       zero in which case\code{NA} is returned.}
\author{Chris Brien}
\examples{
data(exampleData)
Area.WUE <- with(longi.dat, WUI(Area.AGR, Water.Loss))
}
\keyword{data}
\keyword{manip}