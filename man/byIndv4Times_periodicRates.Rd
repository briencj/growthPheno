\name{byIndv4Times_periodicRates}
\alias{byIndv4Times_periodicRates}
\title{Adds the necessary \code{times} and rates values to  a 
             \code{data.frame} to make a set of periodic, or equally spaced, rates.
             }
\description{Adds the times and rates values to a \code{data.frame} to produce 
             a set of periodic rates i.e. rates for equally spaced \code{times}, the 
             difference between  consecutive times being equal to 
             \code{reqd.times.diff}. This assumes that the values for the consecutive 
             differences between the values in the \code{times} column in \code{data} 
             are an integer multiple of the value of \code{reqd.times.diff} and 
             involves inserting times when the integer multiple of a times value is 
             greater than one. It is also assumed that it is appropriate to use the 
             rate value for the time immediately following inserted time(s) for the 
             rate value(s) for the inserted \code{times}. Additional columns can be 
             specified to have values duplicated for the inserted times. 
             
             If \code{avail.time.diffs} is \code{FALSE}, the differences between 
             consecutive time values in \code{data} are calculated.}
\usage{
byIndv4Times_periodicRates(data, responses, 
                           individuals = "Snapshot.ID.Tag", 
                           times = "DAP", 
                           columns2duplicate = NULL, 
                           reqd.times.diff = 1, 
                           avail.times.diffs = FALSE)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns for 
             \code{responses}, \code{individuals} and \code{times}.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} that contains rates that are to be made periodic by inserting   
             extra \code{times} and \code{responses} values to create the data 
             for a set of of equally spaced times separated by 
             \code{reqd.times.diff}.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define(s) the subsets of \code{response} 
             that correspond to the \code{response} values for an individual 
             (e.g. plant, pot, cart, plot or unit). If the column(s) corresponding 
             to \code{individuals} are not \code{\link{factor}}(s) then they will 
             be coerced to \code{\link{factor}}(s). The subsets are formed using 
             \code{\link{split}}.}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times, either as a \code{\link{numeric}}, 
             \code{\link{factor}}, or \code{\link{character}}, at which the data 
             was collected. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
 \item{columns2duplicate}{A \code{\link{character}} giving the name(s) of 
             columns in \code{data} whose values are to be duplicated alongside 
             the inserted \code{times} and \code{responses}. Generally, 
             these should be columns whose values do not differ between the 
             \code{times} values.}
 \item{reqd.times.diff}{A single-valued \code{\link{numeric}} giving the time 
             difference required for the set of equally spaced \code{times}. 
             A natural value for this is the unit of time for the rates specified 
             in \code{responses}. For example, for a rate that measures the 
             change per day, the unit of time is one day and so 1 is the natural 
             value for \code{reqd.times.diff}. All time differences in 
             \code{data} must be an integer multiple of \code{reqd.times.diff} 
             for a set of equally spaced times to be possible.}
 \item{avail.times.diffs}{A \code{\link{logical}} indicating whether there is an 
             appropriate column of \code{times} differences in \code{data} that 
             can be used as the denominator in computing the rates. If \code{TRUE}, 
             it will be assumed that the name of the column is the value of 
             \code{times} with \code{.diffs} appended. If \code{FALSE}, a column, 
             whose column name will be the value of \code{times} with \code{.diffs} 
             appended, will be formed and saved in the result, overwriting any 
             existing columns with the constructed name in \code{data}. It will 
             be calculated using the values of times in \code{data}.}
}
\value{A \code{\link{data.frame}} obtained by expanding \code{data} to include all 
       combinations of the values of \code{individuals} and the set of equally 
       spaced \code{times} between the minimum and maximun values of \code{times},
       the spacing being the value of \code{reqd.times.diff}. 
       For the rate(s) specified in the argument \code{responses}, values of the 
       rate for the \code{times} immediately above an inserted \code{times} value 
       are assigned to the values of the \code{responses} for the inserted 
       \code{times}.}
\author{Chris Brien}
\seealso{\code{\link{byIndv4Times_WaterUse}}}
\examples{
data(exampleData)
longi.dat <- within(longi.dat, WUR <- WU/DAP.diffs)
longi.dat <- byIndv4Times_periodicRates(data = longi.dat, 
                                        responses = "WUR", 
                                        individuals = "Snapshot.ID.Tag", 
                                        times = "DAP", 
                                        avail.times.diffs = TRUE)
}
\keyword{data}
\keyword{manip}