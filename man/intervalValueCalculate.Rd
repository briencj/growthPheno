\name{intervalValueCalculate}
\alias{intervalValueCalculate}
\title{Calculates a single value that is a function of an individual's values for a response over a specified time interval}
\description{Splits the values of a response into subsets corresponding 
             individuals and applies a function that calculates a single 
             value from each individual's observations during a specified 
             time interval. It includes the 
             ability to calculate the observation that corresponds to 
             the calculated value of the function.}
\usage{
intervalValueCalculate(response, weights=NULL, individuals = "Snapshot.ID.Tag",
                       FUN = "max", which.obs = FALSE, which.levels = NULL, 
                       start.time=NULL, end.time=NULL, times.factor = "Days", 
                       suffix.interval=NULL, data, sep=".", na.rm=TRUE, ...)
}
\arguments{
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} from which the values of \code{FUN} are to be  calculated.}
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
 \item{start.time}{A \code{\link{numeric}} giving the times, 
                   in terms of levels of \code{times.factor}, that will give a 
                   single value for each \code{Snapshot.ID.Tag} and that will 
                   be taken as the observation at the start of the interval for 
                   which the growth rate is to be calculated. If 
                   \code{start.time} is \code{NULL}, the interval will start 
                   with the first observation.}
 \item{end.time}{A \code{\link{numeric}} giving the times,
                 in terms of levels of \code{times.factor}, that will give a 
                 single value for each \code{Snapshot.ID.Tag} and that will be 
                 taken as the observation at the end of the interval for which 
                 the growth rate is to be calculated. If \code{end.time} is 
                 \code{NULL}, the interval will end with the last observation.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels should be numeric values stored as 
             characters.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
               appended to \code{response} to form the 
               name of the column containing the calculated values. If it is
               \code{NULL} then nothing will be appended.}
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
       \code{individuals}, containing a column for the \code{individuals}, 
       a column with the values of the function for the \code{individuals},
       and a column with the values of the \code{times.factor}. The name of 
       the column with the values of the function will be result of 
       concatenating the \code{response}, \code{FUN} and, if it is not 
       \code{NULL},  \code{suffix.interval}, each separated by a full stop.}
\author{Chris Brien}
\seealso{\code{\link{intervalGRaverage}}, \code{\link{intervalGRdiff}}, \code{\link{intervalWUI}}, \code{\link{splitValueCalculate}}, \code{\link{getTimesSubset}}}
\examples{
data(exampleData)
Area.smooth.max <- intervalValueCalculate("Area.smooth", 
                                          start.time = 31, end.time = 35, 
                                          suffix.interval = "31to35",
                                          data = longi.dat)
}
\keyword{data}
\keyword{manip}