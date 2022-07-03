\name{byIndv4Intvl_ValueCalc}
\alias{byIndv4Intvl_ValueCalc}
\title{Calculates a single value that is a \code{function} of the values of 
    an individual for a response in a \code{data.frame} in long format 
    over a specified time interval.}
\description{Splits the values of a response into subsets corresponding 
             individuals and applies a function that calculates a single 
             value from each individual's observations during a specified 
             time interval.  It includes the ability to calculate the 
             observation number that is closest to the calculated value of 
             the function and the assocated values of a \code{\link{factor}} or 
             \code{numeric}.}
\usage{
byIndv4Intvl_ValueCalc(data, response, 
                       individuals = "Snapshot.ID.Tag", times = "DAP", 
                       FUN = "max", which.obs = FALSE, which.values = NULL, 
                       addFUN2name = TRUE, 
                       start.time=NULL, end.time=NULL, suffix.interval=NULL, 
                       sep=".", weights=NULL, na.rm=TRUE, ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             function is to be calculated.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} from which the values of \code{FUN} are to be  calculated.}
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
 \item{FUN}{A \code{\link{character}} giving the name of the \code{function} 
             that calculates the value for each subset.}
 \item{which.obs}{A \code{\link{logical}} indicating whether or not to determine 
             the observation number corresponding to the observed value that is 
             closest to the value of the function, in addition to the value of 
             the function itself. That is, \code{FUN} need not return an observed 
             value of the \code{reponse}, e.g. \code{quantile}.}
 \item{which.values}{A \code{\link{character}} giving the name of the 
             \code{\link{factor}} or \code{\link{numeric}} whose values are 
             associated with the \code{response} values and whose value is to be 
             returned for the observation number whose \code{response} value 
             corresponds to the observed value closest to the value of the 
             function. That is, \code{FUN} need not return an observed value of 
             the \code{reponse}, e.g. \code{quantile}. In the case of multiple 
             observed response values satisfying this condition, the value of 
             the \code{which.values} vector for the first of these is 
             returned.}
 \item{addFUN2name}{A \code{\link{logical}} that, if \code{TRUE}, indicates that 
             the \code{FUN} name is to be added to the names of the columns in 
             the \code{\link{data.frame}} returned by \code{byIndv4Intvl_ValueCalc}.}
 \item{start.time}{A \code{\link{numeric}} giving the times, 
             in terms of levels of \code{times.factor}, that will give a 
             single value for each \code{Snapshot.ID.Tag} and that will 
             be taken as the observation at the start of the interval for 
             which a value is to be calculated. If 
             \code{start.time} is \code{NULL}, the interval will start 
             with the first observation. In the case of multiple 
             observed response values satisfying this condition, the first is 
             returned.}
 \item{end.time}{A \code{\link{numeric}} giving the times,
             in terms of levels of \code{times.factor}, that will give a 
             single value for each \code{Snapshot.ID.Tag} and that will be 
             taken as the observation at the end of the interval for which 
             a value is to be calculated. If \code{end.time} is 
             \code{NULL}, the interval will end with the last observation.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
             appended to \code{response} to form the 
             name of the column containing the calculated values. If it is
             \code{NULL} then nothing will be appended.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
             levels of \code{individuals} are combined. This is needed to avoid 
             using a \code{\link{character}} that occurs in a \code{\link{factor}} to delimit 
             levels when the levels of \code{individuals} are combined to 
             identify subsets.}
 \item{weights}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the weights to be supplied as \code{w} 
             to \code{FUN}.}
 \item{na.rm}{A \code{\link{logical}}  indicating whether \code{NA} values 
             should be stripped before the calculation proceeds.}
 \item{... }{allows for arguments to be passed to \code{FUN}.}
}
\value{A \code{\link{data.frame}}, with the same number of rows as there are 
       \code{individuals}, containing a column for the \code{individuals} and 
       a column with the values of the function for the \code{individuals}. 
       It is also possible to determine observaton numbers or the values of 
       another column in \code{data} for the \code{response} values that are 
       closest to the \code{FUN} results, using either or both of 
       \code{which.obs} and \code{which.values}. If \code{which.obs} is 
       \code{TRUE}, a column with observation numbers is included in the  
       \code{\link{data.frame}}. 
       If \code{which.values} is set to the name of a \code{\link{factor}} or 
       a \code{\link{numeric}}, a column containing the levels of that 
       \code{\link{factor}} or the values of that \code{\link{numeric}} is 
       included in the  \code{\link{data.frame}}.
       
       The name of the column with the values of the function will be result of 
       concatenating the \code{response}, \code{FUN} and, if it is not 
       \code{NULL}, \code{suffix.interval}, each separated by a full stop. 
       If \code{which.obs} is \code{TRUE}, the column name for the obervations 
       numbers will have \code{.obs} added after \code{FUN} into the column name 
       for the function values; if \code{which.values} is specified, 
       the column name for these values will have a full stop followed by 
       \code{which.values} added after \code{FUN} into the column name 
       for the function values.
}
\author{Chris Brien}
\seealso{\code{\link{byIndv4Intvl_GRsAvg}}, \code{\link{byIndv4Intvl_GRsDiff}}, 
         \code{\link{byIndv4Intvl_WaterUse}}, \code{\link{splitValueCalculate}}, \cr 
         \code{\link{getTimesSubset}}}
\examples{
data(exampleData)
sPSA.max <- byIndv4Intvl_ValueCalc(data = longi.dat, 
                                   response = "sPSA", times = "DAP", 
                                   start.time = 31, end.time = 35, 
                                   suffix.interval = "31to35")
AGR.max.dat <- byIndv4Intvl_ValueCalc(data = longi.dat, 
                                      response = "sPSA", times = "DAP", 
                                      FUN="max", 
                                      start.time = 31, end.time = 35, 
                                      suffix.interval = "31to35",
                                      which.values = "DAP", 
                                      which.obs = TRUE)
}
\keyword{data}
\keyword{manip}