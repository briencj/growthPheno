\name{intervalGRdiff}
\alias{intervalGRdiff}
\title{Calculates the growth rates for a specified time interval}
\description{Using the values of the responses, calculates the specified 
             combination of the Absolute Growth Rates 
             using differences (AGR), the Proportionate Growth Rates (PGR) 
             and Relative Growth Rates using log differences (RGR) between 
             two nominated time points.

             \strong{Note: this function is soft deprecated and may be removed in 
             future versions. \cr Use \code{\link{byIndv4Intvl_GRsDiff}}.}
}
\usage{
intervalGRdiff(responses, individuals = "Snapshot.ID.Tag", 
               which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
               times = "Days", start.time, end.time, suffix.interval, 
               data)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the growth rates are to be  calculated.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (e.g. plant, pot, cart, plot or unit).}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
            to be calculated. It should be a combination of one or more of 
            \code{"AGR"}, \code{"PGR"} and \code{"RGR"}.}
 \item{suffices.rates}{A \code{\link{character}} giving the characters to be 
                    appended to the names of the responses in constructing the 
                    names of the columns containing the calculated growth rates.            
                    The order of the suffices in \code{suffices.rates} should 
                    correspond to the order of the elements of 
                    \code{which.rates}.}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used in calculating growth rates and,  
             if a \code{\link{factor}} or \code{\link{character}}, the values should 
             be numerics stored as characters.}
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
\code{response} at the \code{end.time}  and \code{start.time} divided by the 
difference between \code{end.time}  and \code{start.time}. 
The PGR is calculated as the ratio of \code{response} at the \code{end.time} 
to that at \code{start.time} and the ratio raised to the power of the 
reciprocal of the difference between \code{end.time}  and \code{start.time}. 
The \code{RGR} is calculated as the \code{log} of the PGR and so is equal to 
the difference between the logarithms of \code{response} 
at the \code{end.time}  and \code{start.time} divided by the difference 
between \code{end.time}  and \code{start.time}.
}
\author{Chris Brien}
\seealso{\code{\link{intervalGRaverage}}, \code{\link{intervalWUI}},  \code{\link{getTimesSubset}}, \code{\link{GrowthRates}}, \cr
         \code{\link{splitSplines}}, \code{\link{splitContGRdiff}}}
\examples{
data(exampleData)
sPSA.GR <- intervalGRdiff(responses = "sPSA", times = "DAP", 
                          which.rates = c("AGR","RGR"), 
                          start.time = 31, end.time = 35,
                          suffix.interval = "31to35",
                          data = longi.dat)
}
\keyword{data}
\keyword{manip}