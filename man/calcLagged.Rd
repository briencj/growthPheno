\name{calcLagged}
\alias{calcLagged}
\title{Replaces the values in a vector with the result of applying an operation to it 
      and a lagged value}
\description{Replaces the values in \code{x} with the result of applying an 
             \code{operation} to it and the value that is \code{lag} positions 
             either before it or after it in \code{x}, depending on whether 
             \code{lag} is positive or negative. For positive \code{lag} 
             the first \code{lag} values will be \code{NA}, while for negative 
             \code{lag} the last \code{lag} values will be \code{NA}. 
             When \code{operation} is \code{NULL}, the values are moved 
                  \code{lag} positions down the vector.}
\usage{calcLagged(x, operation = NULL, lag = 1)
}
\arguments{
 \item{x}{A \code{\link{vector}} containing the values on which the calculations are to be made.}
 \item{operation}{A \code{\link{character}} giving the operation to be performed on pairs of 
                  values in \code{x}. If \code{operation} is \code{NULL} then the values are moved 
                  \code{lag} positions down the vector.}
 \item{lag}{A \code{integer} specifying, for the second value in the pair to be operated on, 
            the number positions it is ahead of or behind the current value.}
} 
\value{A \code{\link{vector}} containing  the result of applying \code{operation} to 
       values in \code{x}. For positive \code{lag} the first \code{lag} values will 
       be \code{NA}, while for negative \code{lag} the last \code{lag} values will be \code{NA}.}
\author{Chris Brien}
\seealso{\code{\link{Ops}}}
\examples{
data(exampleData)
longi.dat$DAP.diffs <- calcLagged(longi.dat$xDAP, operation ="-")
}
\keyword{data}
\keyword{manip}