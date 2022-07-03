\name{intervalWUI}
\alias{intervalWUI}
\title{Calculates water use indices (WUI) over a specified time interval to a 
       data.frame}
\description{Calculates the Water Use Index (WUI) between two time points 
             for a set of responses.

             \strong{Note: this function is soft deprecated and may be removed in 
             future versions. \cr Use \code{\link{byIndv4Intvl_WaterUse}}.}
}
\usage{
intervalWUI(responses, water.use = "Water.Use", 
            individuals = "Snapshot.ID.Tag", times = "Days", 
            start.time, end.time, suffix.interval = NULL, 
            data, include.total.water = FALSE, na.rm = FALSE)}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the growth rates are to be calculated.}
 \item{water.use}{A \code{\link{character}} giving the names of the column in 
             \code{data} which contains the water use values.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (plant/cart/plot/unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used identifying the intervals and,  
             if a \code{\link{factor}} or \code{\link{character}}, the values should 
             be numerics stored as characters.}
 \item{start.time}{A \code{\link{numeric}} giving the times,
                   in terms of values in \code{times}, that will give a 
                   single value for each \code{Snapshot.ID.Tag} and that will 
                   be taken as the observation at the start of the interval for 
                   which the WUI is to be calculated.}
 \item{end.time}{A \code{\link{numeric}} giving the times, 
                 in terms of values \code{times}, that will give a 
                 single value for each \code{Snapshot.ID.Tag} and that will be 
                 taken as the observation at the end of the interval for which 
                 the WUI is to be calculated.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
               appended to \code{response} to form the names of the columns 
               containing the calculated the growth rates.}
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             growth rates are to be calculated.}
 \item{include.total.water}{A \code{\link{logical}} indicating whether or not to 
             include a column in the results for the total of \code{water.use} 
             for the interval for each individual.}
 \item{na.rm}{A \code{\link{logical}}  indicating whether \code{NA} values 
              should be stripped before the calculation proceeds.}
}
\value{A \code{\link{data.frame}} containing the WUIs, the name of each 
       column being the concatenation of one of \code{responses},  
       \code{WUI} and, if not \code{NULL}, \code{suffix.interval}, the three 
       components being separated by a full stop. If the total water is to be 
       included, the name of the column will be the concatenation of 
       \code{water.use}, \code{Total} and the suffix, each separated by a full
       stop(`.').}
\details{
The WUI is calculated as the difference between the values of a \code{response} 
at the \code{end.time}  and \code{start.time} divided by the sum of the water use  
after \code{start.time} until \code{end.time}. Thus, the water use up to
\code{start.time} is not included.
}
\author{Chris Brien}
\seealso{\code{\link{intervalGRaverage}}, \code{\link{intervalGRdiff}}, \code{\link{splitValueCalculate}}, \code{\link{getTimesSubset}}, 
         \code{\link{GrowthRates}}}
\examples{
data(exampleData)
PSA.WUI <- intervalWUI(response = "PSA", water.use = "WU", 
                       times = "DAP", start.time = 31, end.time = 35, 
                       suffix = "31to35", 
                       data = longi.dat, include.total.water = TRUE)
}
\keyword{data}
\keyword{manip}