\name{getTimesSubset}
\alias{getTimesSubset}
\title{Forms a subset of \code{responses} in \code{data} that contains their values for the nominated times}
\description{Forms a subset of \code{responses} in \code{data} that contains their values for the 
             nominated times.}
\usage{
getTimesSubset(responses, times.factor = "Days", data, which.times, 
               suffix = NULL, include.times.factor = FALSE, 
               include.individuals = FALSE, individuals = "Snapshot.ID.Tag")}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} whose values are to be subsetted.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used to identify the subset and 
             should be numeric values stored as characters.}
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             growth rates are to be calculated.}
 \item{which.times}{A \code{\link{vector}} giving the times that are to be selected.}
 \item{suffix}{A \code{\link{character}} giving the suffix to be appended to 
               \code{responses} to form the names of the columns containing the 
               subset.}
 \item{include.times.factor}{A \code{\link{logical}} indicating whether or not 
             to include the \code{times.factor} in the result, the name in the result 
             having the \code{suffix} with a separating full appended.}
 \item{include.individuals}{A \code{\link{logical}} indicating whether or not 
             to include the \code{individuals} column in the result.}
 \item{individuals}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing an identifier for each individual.}
}
\value{A \code{\link{data.frame}} containing the subset of \code{responses} ordered 
       by as many of the initial columns as are required to uniquely identify each row 
       (see \code{\link{order}} for more information). The names of the columns for 
       \code{responses} and \code{times.factor} in the subset are the concatenation of   
       their names in \code{data} and \code{suffix} separated by a full stop.}
\author{Chris Brien}
\examples{
data(exampleData)
AreaLast <- getTimesSubset("Area.smooth", data = longi.dat, 
                           which.times = c(42), suffix = "last")
}
\keyword{data}
\keyword{manip}