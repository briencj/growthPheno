\name{intervalGRdiff}
\alias{intervalGRdiff}
\title{Calculates the growth rates for a specified time interval}
\description{Using the values of the responses, calculates the specified 
             combination of the Absolute Growth Rates 
             using differences (AGR), the Proportionate Growth Rates (PGR) 
             and Relative Growth Rates using log differences (RGR) between 
             two nominated time points.}
\usage{
intervalGRdiff(responses, individuals = "Snapshot.ID.Tag", 
               which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
               times.factor = "Days", start.times, end.times, suffix.interval, 
               data)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the growth rates are to be  calculated.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
                    \code{factor}(s) that define the subsets of the \code{data} 
                    for which each subset corresponds to the \code{responses} 
                    for an individual.}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
                    to be calculated. It should be a combination 
                    \code{"AGR"}, \code{"PGR"} and \code{"RGR"}.}
 \item{suffices.rates}{A \code{\link{character}} giving the characters to be 
                    appended to the names of the responses in constructing the 
                    names of the columns containing the calculated growth rates.            
                    The order of the suffices in \code{suffices.rates} should 
                    correspond to the order of the elements of 
                    \code{which.rates}.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used in calculating growth rates and 
             should be numeric values stored as characters.}
 \item{start.times}{A \code{\link{numeric}} giving the times,
                    in terms of levels of \code{times.factor}, that will give a 
                    single value for each \code{Snapshot.ID.Tag} and that will 
                    be taken as the observation at the start of the interval for 
                    which the growth rate is to be calculated.}
 \item{end.times}{A \code{\link{numeric}} giving the times, 
                  in terms of levels of \code{times.factor}, that will give a 
                  single value for each \code{Snapshot.ID.Tag} and that will be 
                  taken as the observation at the end of the interval for which 
                  the growth rate is to be calculated.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
               appended to \code{response} to form the names of the columns 
               containing the calculated the growth rates.}
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             growth rates are to be calculated.}
}
\value{A \code{\link{data.frame}} with the growth rates. 
       The name of each column is the concatenation of (i) one of 
       \code{responses}, (ii) one of \code{AGR}, \code{PGR} or \code{RGR},
       or the appropriate element of \code{suffices.rates}, and (iii) 
       \code{suffix.interval}, the three components being separated by 
        full stops.}
\details{
The \code{AGR} is calculated as the difference between the values of 
\code{response} at the \code{end.times}  and \code{start.times} divided by the 
difference between \code{end.times}  and \code{start.times}. 
The PGR is calculated as the ratio of \code{response} at the \code{end.times} 
to that at \code{start.times} and the ratio raised to the power of the 
reciprocal of the difference between \code{end.times}  and \code{start.times}. 
The \code{RGR} is calculated as the \code{log} of the PGR and so is equal to 
the difference between the logarithms of \code{response} 
at the \code{end.times}  and \code{start.times} divided by the difference 
between \code{end.times}  and \code{start.times}.
}
\author{Chris Brien}
\seealso{\code{\link{intervalGRaverage}}, \code{\link{intervalWUI}},  \code{\link{getTimesSubset}}, \code{\link{GrowthRates}}, \code{\link{splitSplines}}, \code{\link{splitContGRdiff}}}
\examples{
data(exampleData)
Area.smooth.GR <- intervalGRdiff("Area.smooth", which.rates = c("AGR","RGR"), 
                                 start.times = 31, end.times = 35,
                                 suffix.interval = "31to35",
                                 data = longi.dat)
}
\keyword{data}
\keyword{manip}