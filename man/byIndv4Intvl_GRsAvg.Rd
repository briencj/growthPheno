\name{byIndv4Intvl_GRsAvg}
\alias{byIndv4Intvl_GRsAvg}
\title{Calculates the growth rates for a specified time interval for individuals 
       in a \code{data.frame} in long format by taking weighted averages of 
       growth rates for times within the interval.}
\description{Using previously calculated growth rates over time, calculates 
       the Absolute Growth Rates for a specified interval using the weighted 
       averages of AGRs for each time point in the interval (AGR) and the 
       Relative Growth Rates for a specified interval using the weighted 
       geometric means of RGRs for each time point in the interval (RGR).}
\usage{
byIndv4Intvl_GRsAvg(data, responses, 
                    individuals = "Snapshot.ID.Tag", times = "DAP", 
                    which.rates = c("AGR","RGR"), 
                    suffices.rates=c("AGR","RGR"), 
                    start.time, end.time, suffix.interval, 
                    sep=".", na.rm=TRUE)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns from which the 
             growth rates are to be calculated.}
 \item{responses}{A \code{\link{character}} giving the names of the 
                  \code{responses} for which there are columns 
                  in \code{data} that contain the growth rates that are to be 
                  averaged. The names of the growth rates should have 
                  either \code{AGR} or \code{RGR} appended to the 
                  \code{responses} names.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (plant/cart/plot/unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used in calculating growth rates and,  
             if a \code{\link{factor}} or \code{\link{character}}, the values should 
             be numerics stored as characters.}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
             to be averaged to obtain growth rates for an interval. 
             It should be a combination of one or more of 
             \code{"AGR"} and \code{"RGR"}.}
 \item{suffices.rates}{A \code{\link{character}} giving the suffices to be 
               appended to \code{response} to form the names of the columns 
               containing the calculated the growth rates and in which growth 
               rates are to be stored. Their elements will be matched with 
               those of \code{which.rates}.}
 \item{start.time}{A \code{\link{numeric}} giving the times,
                   in terms of values in \code{times}, that will give a 
                   single value for each \code{Snapshot.ID.Tag} and that will 
                   be taken as the observation at the start of the interval for 
                   which the growth rate is to be calculated.}
 \item{end.time}{A \code{\link{numeric}} giving the times, 
                 in terms of values \code{times}, that will give a 
                 single value for each \code{Snapshot.ID.Tag} and that will be 
                 taken as the observation at the end of the interval for which 
                 the growth rate is to be calculated.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
               appended to \code{response}.\code{suffices.rates} to form the 
               names of the columns containing the calculated the growth rates.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a \code{\link{factor}} to delimit 
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
obtained from the \code{times}. They are taken as the sum of half the time subintervals before and after each time, except for the end points; the end points are taken to be the subintervals at the start and end of the interval.
}
\author{Chris Brien}
\seealso{\code{\link{byIndv4Intvl_GRsDiff}}, \code{\link{byIndv4Intvl_WaterUse}}, \code{\link{splitValueCalculate}}, \code{\link{getTimesSubset}}, \cr 
         \code{\link{GrowthRates}}, \code{\link{byIndv4Times_SplinesGRs}}, \code{\link{splitContGRdiff}}}
\examples{
data(exampleData)
longi.dat <- byIndv4Times_SplinesGRs(data = longi.dat, 
                                     response="PSA", response.smoothed = "sPSA", 
                                     individuals = "Snapshot.ID.Tag", 
                                     times = "DAP", 
                                     df = 4, 
                                     rates.method = "deriv", 
                                     which.rates = c("AGR", "RGR"), 
                                     suffices.rates = c("AGRdv", "RGRdv"))
sPSA.GR <- byIndv4Intvl_GRsAvg(data = longi.dat, 
                               response="sPSA", times = "DAP", 
                               which.rates = c("AGR","RGR"), 
                               suffices.rates = c("AGRdv","RGRdv"), 
                               start.time = 31, end.time = 35, 
                               suffix.interval = "31to35")
}
\keyword{data}
\keyword{manip}