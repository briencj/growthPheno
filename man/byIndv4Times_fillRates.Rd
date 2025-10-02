\name{byIndv4Times_fillRates}
\alias{byIndv4Times_fillRates}
\title{Adds to the times and rates values in  a \code{data.frame} to produce 
       a set of equally spaced times with values for the rates  for the inserted times.}
\description{Adds to the values in a \code{data.frame} for the \code{times} column 
             and the rates nominated in \code{responses} to produce a set of 
             equally spaced times with values for the rates for the inserted times. 
             
             If \code{avail.time.diffs} is \code{FALSE}, the differences between 
             consecutive time values are  calculated.}
\usage{
byIndv4Times_fillRates(data, responses, 
                       individuals = "Snapshot.ID.Tag", times = "DAP", 
                       avail.times.diffs = FALSE, 
                       ntimes2span = 2)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns for 
             \code{responses}, \code{individuals} and \code{times}.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for which rates are to be inserted at the values of  
             \code{times} that have been added so that there are values for the set 
             of equally spaced times.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of \code{response} 
             that correspond to the \code{response} values for an individual 
             (e.g. plant, pot, cart, plot or unit) for which growth rates are 
             to be calculated continuously. If the columns corresponding to 
             \code{individuals} are not \code{\link{factor}}(s) then they will 
             be coerced to \code{\link{factor}}(s). The subsets are formed using 
             \code{\link{split}}.}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, 
             or \code{\link{character}}. It will be used in calculating the 
             growth rates. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
 \item{avail.times.diffs}{A \code{\link{logical}} indicating whether there is an 
             appropriate column of \code{times} differences that can be used as 
             the denominator in computing the growth rates. If \code{TRUE}, it will 
             be assumed that the name of the column is the value of  \code{times} 
             with \code{.diffs} appended. If \code{FALSE}, a column, whose 
             column name will be the value of \code{times} with \code{.diffs} 
             appended, will be formed and saved in the result, overwriting any 
             existing columns with the constructed name in \code{data}. It will 
             be calculated using the values of times in \code{data}.}
 \item{ntimes2span}{A \code{\link{numeric}} giving the number of values in 
             \code{times} to span in calculating growth rates by differencing. 
             Each growth rate is calculated as the difference in the values of 
             one of the \code{responses} for pairs of \code{times} values that 
             are spanned by \code{ntimes2span} \code{times} values divided by 
             the difference between this pair of \code{times} values. For 
             \code{ntimes2span} set to 2, a growth rate is the  
             difference between consecutive pairs of values of one of the 
             \code{responses} divided by the difference between consecutive 
             pairs of \code{times} values.}
}
\value{If the differences between consecutive values of \code{times} in \code{data} 
       are not equal, the \code{\link{data.frame}} is expanded to include all 
       combinations of the values of \code{individuals} and  the set of equally 
       spaced \code{times} between the minimum and maximun values of \code{times}. 
       For the rate(s) specified in the argument \code{responses}, values of the 
       rate for the \code{times} immediately above an inserted \code{times} value 
       are assigned to the values of the \code{responses} for the inserted 
       \code{times}.}
\author{Chris Brien}
\seealso{\code{\link{byIndv4Times_WaterUse}}}
\examples{
data(exampleData)
longi.dat <- within(longi.dat, WUR <- WU/DAP.diffs)
longi.dat <- byIndv4Times_fillRates(data = longi.dat, 
                                    response = "WUR", 
                                    individuals = "Snapshot.ID.Tag", 
                                    times = "DAP", 
                                    avail.times.diffs = TRUE)
}
\keyword{data}
\keyword{manip}