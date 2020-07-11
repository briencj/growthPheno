\name{intervalGRaverage}
\alias{intervalGRaverage}
\title{Calculates the growth rates for a specified time interval by taking weighted averages of growth rates for times within the interval}
\description{Using previously calculated growth rates over time, calculates 
             the Absolute Growth Rates for a specified interval 
             using the weighted averages of AGRs for each time point in the 
             interval (AGR) and the Relative Growth Rates for a specified 
             interval using the weighted geometric means of RGRs for each time 
             point in the interval (RGR).}
\usage{
intervalGRaverage(responses, individuals = "Snapshot.ID.Tag", 
                  which.rates = c("AGR","RGR"), suffices.rates=c("AGR","RGR"), 
                  start.time, end.time, times.factor = "Days", suffix.interval, 
                  data, sep=".", na.rm=TRUE)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the 
                  \code{responses} for which there are columns 
                  in \code{data} that contain the growth rates that are to be 
                  averaged. The names of the growth rates should have 
                  either \code{AGR} or \code{RGR} appended to the 
                  \code{responses} names.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
                    \code{factor}(s) that define the subsets of the \code{data} 
                    for which each subset corresponds to the \code{responses} 
                    for an individual.}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
                    to be averaged to obtain growth rates for an interval. 
                    It should be a combination \code{"AGR"} and \code{"RGR"}.}
 \item{suffices.rates}{A \code{\link{character}} giving the suffices to be 
               appended to \code{response} to form the names of the columns 
               containing the calculated the growth rates and in which growth 
               rates are to be stored. Their elements will be matched with 
               those of \code{which.rates}.}
 \item{start.time}{A \code{\link{numeric}} giving the times,
                   in terms of levels of \code{times.factor}, that will give a 
                   single value for each \code{Snapshot.ID.Tag} and that will 
                   be taken as the observation at the start of the interval for 
                   which the growth rate is to be calculated.}
 \item{end.time}{A \code{\link{numeric}} giving the times, 
                 in terms of levels of \code{times.factor}, that will give a 
                 single value for each \code{Snapshot.ID.Tag} and that will be 
                 taken as the observation at the end of the interval for which 
                 the growth rate is to be calculated.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used in calculating growth rates and 
             should be numeric values stored as characters.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
               appended to \code{response}.\code{suffices.rates} to form the 
               names of the columns containing the calculated the growth rates.}
 \item{data}{A \code{\link{data.frame}} containing the columns from which the 
             growth rates are to be calculated.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a factor to delimit 
            levels when the levels of \code{individuals} are combined to 
            identify subsets.}
 \item{na.rm}{A \code{\link{logical}} indicating whether \code{NA} values should be 
              stripped before the calculation of weighted means proceeds.}
}
\value{A \code{\link{data.frame}} with the growth rates. 
       The name of each column is the concatenation of (i) one of 
       \code{responses}, (ii) one of \code{AGR}, \code{PGR} or \code{RGR},
       or the appropriate element of \code{suffices.rates}, and (iii) 
       \code{suffix.interval}, the three components being separated by 
        full stops.}
\details{
The \code{AGR} for an interval is calculated as the weighted mean of the 
AGRs for times within the interval. The RGR is calculated as the weighted 
geometric mean of the RGRs for times within the interval; in fact the exponential is taken of the weighted means of the logs of the RGRs. The weights are 
obtained from the \code{times.factor}. They are taken as the sum of half the time subintervals before and after each time, except for the end points; the end points are taken to be the subintervals at the start and end of the interval.
}
\author{Chris Brien}
\seealso{\code{\link{intervalGRdiff}}, \code{\link{intervalWUI}}, \code{\link{splitValueCalculate}}, \code{\link{getTimesSubset}}, \code{\link{GrowthRates}}, \code{\link{splitSplines}}, \code{\link{splitContGRdiff}}}
\examples{
data(exampleData)
longi.dat <- splitSplines(longi.dat, response="Area", x="xDays", 
                          INDICES = "Snapshot.ID.Tag", 
                          df = 4, deriv=1, suffices.deriv="AGRdv", RGR="RGRdv")
Area.smooth.GR <- intervalGRaverage("Area.smooth", which.rates = c("AGR","RGR"), 
                                    suffices.rates = c("AGRdv","RGRdv"), 
                                    start.time = 31, end.time = 35, 
                                    suffix.interval = "31to35",
                                    data = longi.dat)
}
\keyword{data}
\keyword{manip}