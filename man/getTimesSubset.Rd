\name{getTimesSubset}
\alias{getTimesSubset}
\title{Forms a subset of \code{responses} in \code{data} that contains their values for the nominated times}
\description{Forms a subset of each of the \code{responses} in \code{data} that contains their values for the 
             nominated times in a single column.}
\usage{
getTimesSubset(data, responses, 
               individuals = "Snapshot.ID.Tag", times = "DAP", 
               which.times, suffix = NULL, sep.suffix.times = ".", 
               include.times = FALSE, include.individuals = FALSE)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             growth rates are to be calculated.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} whose values are to be subsetted.}
 \item{individuals}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing an identifier for each individual 
             (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used to identify the subset and,  
             if a \code{\link{factor}} or \code{\link{character}}, the values should 
             be numerics stored as characters.}
 \item{which.times}{A \code{\link{vector}} giving the times that are to be selected.}
 \item{suffix}{A \code{\link{character}} giving the suffix to be appended to 
               \code{responses} to form the names of the columns containing the 
               subset.}
 \item{sep.suffix.times}{A \code{\link{character}} giving the separator to use in 
            appending a suffix for times to a trait. For no separator, set to 
            \code{""}.}
 \item{include.times}{A \code{\link{logical}} indicating whether or not 
             to include the \code{times} in the result, the name in the result 
             having the \code{suffix} with a separating full appended.}
 \item{include.individuals}{A \code{\link{logical}} indicating whether or not 
             to include the \code{individuals} column in the result.}
}
\value{A \code{\link{data.frame}} containing the subset of \code{responses} ordered 
       by as many of the initial columns of \code{data} as are required to uniquely 
       identify each row (see \code{\link{order}} for more information). The names 
       of the columns for each of the \code{responses} and for \code{times} in the 
       subset are the concatenation of their names in \code{data} and \code{suffix}, 
       separated by a full stop.}
\author{Chris Brien}
\examples{
data(exampleData)
sPSALast <- getTimesSubset("sPSA", data = longi.dat, times = "DAP", 
                           which.times = c(42), suffix = "last")
}
\keyword{data}
\keyword{manip}