\name{byIndv4Times_GRsDiff}
\alias{byIndv4Times_GRsDiff}
\title{Adds, to a \code{data.frame}, the growth rates calculated for consecutive times 
       for individuals in a \code{data.frame} in long format by differencing 
       response values.}
\description{Uses \code{\link{AGRdiff}}, \code{\link{PGR}} and 
             \code{\link{RGRdiff}} to calculate growth rates continuously 
             over time for the \code{response} by differencing pairs of pairs 
             of \code{response} values  and 
             stores the results in \code{data}. The subsets are those values 
             with the same levels combinations of the \code{\link{factor}}s listed in 
             \code{individuals}.}
\usage{
byIndv4Times_GRsDiff(data, responses, 
                     individuals = "Snapshot.ID.Tag", times = "DAP", 
                     which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                     avail.times.diffs = FALSE, ntimes2span = 2)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns for which 
             growth rates are to be calculated.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for which growth rates are to be calculated.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of \code{response} 
             that correspond to the \code{response} values for an individual 
             (plant/cart/plot/unit) for which growth rates are to be calculated 
             continuously. If the columns corresponding to \code{individuals} 
             are not \code{\link{factor}}(s) then they will be coerced to 
             \code{\link{factor}}(s). The subsets are formed using 
             \code{\link{split}}.}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, 
             or \code{\link{character}}. It will be used in calculating the 
             growth rates. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
             to be calculated. It should be a combination of one or more of 
             \code{"AGR"}, \code{"PGR"} and \code{"RGR"}.}
 \item{suffices.rates}{A \code{\link{character}} giving the characters to be 
             appended to the names of the responses to provide the names of the 
             columns containing the calculated growth rates. The order of the 
             suffices in \code{suffices.rates} should correspond to the order 
             of the elements of \code{which.rates}. If \code{NULL}, the values 
             of \code{which.rates} are used.}
 \item{avail.times.diffs}{A \code{\link{logical}} indicating whether there is an 
             appropriate column of \code{times} diffserences that can be used as 
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
\value{A \code{\link{data.frame}} containing  \code{data} to which has been 
       added i) a column for the differences between the \code{times}, 
       if it is not already in \code{data}, and (ii) columns with growth rates. 
       The name of the column for \code{times} differences will be the value 
       of \code{times} with \code{".diffs"} appended. The name for each of 
       the growth-rate columns will be either the value of \code{response} 
       with one of \code{".AGR"}, \code{".PGR"} or \code{"RGR"}, or the 
       corresponding value from \code{suffices.rates} appended. Each growth 
       rate will be positioned at observation 
       \code{ceiling(ntimes2span + 1) / 2} relative to the two times from 
       which the growth rate is calculated.}
\author{Chris Brien}
\seealso{\code{\link{smoothSpline}}, \code{\link{byIndv4Times_SplinesGRs}}}
\examples{
data(exampleData)
longi.dat <- byIndv4Times_GRsDiff(data = longi.dat, 
                                  response = "sPSA", 
                                  individuals = "Snapshot.ID.Tag", 
                                  times = "DAP", 
                                  which.rates=c("AGR", "RGR"), 
                                  avail.times.diffs = TRUE)
}
\keyword{data}
\keyword{manip}