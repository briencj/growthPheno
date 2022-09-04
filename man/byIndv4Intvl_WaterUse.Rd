\name{byIndv4Intvl_WaterUse}
\alias{byIndv4Intvl_WaterUse}
\title{Calculates, water use traits (WU, WUR, WUI) over a specified time interval 
       for each individual in a \code{data.frame} in long format.}
\description{Calculates one or more of water use (WU), water use rate (WUR), and, 
       for a set of responses,  water use indices (WUI)s over a specified time 
       interval for each individual in a \code{data.frame} in long format.}
\usage{
byIndv4Intvl_WaterUse(data, water.use = "Water.Use", responses = NULL, 
                      individuals = "Snapshot.ID.Tag", times = "DAP", 
                      trait.types = c("WU", "WUR", "WUI"), 
                      suffix.rate = "R", suffix.index = "I", 
                      sep.water.traits = "", sep.responses = ".", 
                      start.time, end.time, 
                      suffix.interval = NULL, sep.suffix.interval = ".", 
                      na.rm = FALSE)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             water use traits are to be calculated.}
 \item{water.use}{A \code{\link{character}} giving the names of the column in 
             \code{data} that contains the water use values.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for which \code{WUI}s are to be calculated.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
             \code{\link{factor}} that defines the subsets of the \code{data} 
             for which each subset corresponds to the \code{response} values for 
             an individual (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used identifying the intervals and,  
             if a \code{\link{factor}} or \code{\link{character}}, the values should 
             be numerics stored as characters.}
 \item{trait.types}{A \code{\link{character}} listing the trait types to compute
             and return. It should be some combination of \code{WU}, \code{WUR} 
             and \code{WUI}, or be \code{all}. See Details for how each 
             is calculated.}
 \item{suffix.rate}{A \code{\link{character}} giving the label to be appended 
             to the value of \code{water.use} to form the name of the WUR.}
 \item{suffix.index}{A \code{\link{character}} giving the label to be appended 
             to the value of \code{water.use} to form the name of the WUI.}
 \item{sep.water.traits}{A \code{\link{character}} giving the character(s) to be 
             used to separate the \code{suffix.rate} and \code{suffix.index} values 
             from the \code{responses} values in constructing the name for a new 
             rate/index. The default of \code{""} results in no separator.}
 \item{sep.responses}{A \code{\link{character}} giving the character(s) to be used 
             to separate the \code{suffix.rate} value from a\code{responses} 
             value in constructing the name for a new index. For no separator, 
             set to \code{""}.}
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
             appended to the names of the columns for the water use traits to 
             indicate the interval for which the  traits have been calculated .}
 \item{sep.suffix.interval}{A \code{\link{character}} giving the separator to use in 
            appending \code{suffix.inteval} to a growth rate. For no separator, set to 
            \code{""}.}
 \item{na.rm}{A \code{\link{logical}}  indicating whether \code{NA} values 
             should be stripped before the calculation proceeds.}
}
\value{A \code{\link{data.frame}} containing the \code{individuals} column, \code{WU} 
       and/or \code{WUR} and, if requested, a \code{WUI} for each 
       element of \code{responses}. The names of \code{WU} and \code{WUR} will have 
       \code{suffix.interval} appended, if it is not \code{NULL}, separated by a 
       full stop (`.'). The name of each \code{WUI} will be the concatenation of an 
       element of \code{responses} with \code{WUI} and, if not \code{NULL}, 
       \code{suffix.interval}, the three components being separated by a 
       full stop (`.').}
\details{
\describe{
\item{\code{WU}}{ is the water use and is the sum of the water use after \code{start.time} 
until \code{end.time}. Thus, the water use up to \code{start.time} is not included.} 
\item{\code{WUR}}{ is the Water Use Rate and is \code{WU} divided by the difference 
between \code{end.time} and \code{start.time}.}
\item{\code{WUI}}{ is the Water Use Index and is calculated as a \code{response} difference 
between the start.time and the end.time, which s then divided by the \code{WU}.}
}}
\author{Chris Brien}
\seealso{\code{\link{byIndv4Intvl_GRsAvg}}, \code{\link{byIndv4Intvl_GRsDiff}}, 
          \code{\link{splitValueCalculate}}, \code{\link{getTimesSubset}}, \cr
         \code{\link{GrowthRates}}}
\examples{
data(exampleData)
WU.WUI_31_35 <- byIndv4Intvl_WaterUse(data = longi.dat, 
                                      water.use = "WU", responses = "PSA", 
                                      times = "DAP", 
                                      trait.types = c("WUR","WUI"), 
                                      suffix.rate = ".Rate", 
                                      suffix.index = ".Index",
                                      start.time = 31, end.time = 35, 
                                      suffix.interval = "31to35")
}
\keyword{data}
\keyword{manip}