\name{splitValueCalculate}
\alias{splitValueCalculate}
\title{Calculates a single value that is a function of an individual's values for a response}
\description{Splits the values of a response into subsets corresponding 
             individuals and applies a function that calculates a single 
             value to each individual's observations.
             It includes the ability to calculate the observation that 
             corresponds to the calculated value of the function.}
\usage{
splitValueCalculate(response, weights=NULL, individuals = "Snapshot.ID.Tag", 
                    FUN = "max", which.obs = FALSE, which.levels = NULL, 
                    data, na.rm=TRUE, sep=".", ...)}
\arguments{
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} from which the values of \code{FUN} are to be  
              calculated.}
 \item{weights}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the weights to be supplied as \code{w} 
             to \code{FUN}.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{factor}(s) that define the subsets of the \code{data} for 
             which each subset corresponds to the \code{response} value for an 
             individual.}
 \item{FUN}{A \code{\link{character}} giving the name of the \code{function} 
             that calculates the value for each subset.}
 \item{which.obs}{A \code{\link{logical}} indicating whether or not to determine 
             the observation corresponding to the value of the function, 
             instead of the value of the function itself.}
 \item{which.levels}{A \code{\link{character}} giving the name of the 
                     \code{\link{factor}} whose levels are to be identified as 
                     the level of the observation whose value matches the 
                     value of the function.}
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             function is to be calculated.}
 \item{na.rm}{A \code{\link{logical}}  indicating whether \code{NA} values 
              should be stripped before the calculation proceeds.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a factor to delimit 
            levels when the levels of \code{individuals} are combined to 
            identify subsets.}
 \item{... }{allows for arguments to be passed to \code{FUN}.}
}
\value{A \code{\link{data.frame}}, with the same number of rows as there are 
       \code{individuals}, containing the values of the function for the 
       \code{individuals}.}
\author{Chris Brien}
\seealso{\code{\link{splitContGRdiff}}, \code{\link{splitSplines}}}
\examples{
data(exampleData)
Area.smooth.max <- splitValueCalculate("Area.smooth", data = longi.dat)
}
\keyword{data}
\keyword{manip}