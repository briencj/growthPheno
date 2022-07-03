\name{cumulate}
\alias{cumulate}
\title{Calculates the cumulative sum, ignoring the first element if exclude.1st is TRUE}
\description{Uses \code{cumsum} to calculate the cumulative sum, ignoring the first element 
             if exclude.1st is \code{TRUE}.}
\usage{cumulate(x, exclude.1st = FALSE, na.rm = FALSE, ...)}
\arguments{
 \item{x}{A \code{\link{vector}} containing the values to be cumulated.}
 \item{exclude.1st}{A \code{\link{logical}} indicating whether or not the first value of 
                    the cumulative sum is to be \code{NA}.}
 \item{na.rm}{A \code{\link{logical}} indicating whether \code{NA} values should be 
              stripped before the computation proceeds}
 \item{...}{allows passing of arguments to other functions; not used at present.}
} 
\value{A \code{\link{vector}} containing  the cumulative sum.}
\author{Chris Brien}
\seealso{\code{\link{cumsum}}}
\examples{
data(exampleData)
PSA.cum <- cumulate(longi.dat$PSA)
}
\keyword{data}
\keyword{manip}