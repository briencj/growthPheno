\name{traitExtractFeatures}
\alias{traitExtractFeatures}
\title{Extract features, that are single-valued for each individual, from traits observed over time.}
\description{Extract one or more sets of features from traits observed over time, the result being 
             traits that have a single value for each individual. The sets of features are:
             \enumerate{
               \item{\bold{single times} -- the value for each individual for a single time. 
                     (uses \code{getTimesSubset})}
               \item{\bold{growth rates for a time interval} -- the average growth rate (AGR 
                     and/or RGR) over a time interval for each individual. 
                     (uses \code{byIndv4Intvl_GRsDiff} or \code{byIndv4Intvl_GRsAvg})}
               \item{\bold{water use traits for a time interval} -- the total water use (WU), 
                     the water use rate (WUR) and the water use index (WUI) over a time 
                     interval for each individual. (uses \code{byIndv4Intvl_WaterUse})}
               \item{\bold{whole of imaging period} -- the total over the whole imaging period 
                     of a trait for each individual. 
                     (uses byIndv4Intvl_ValueCalc)}
               \item{\bold{maximum} -- the maximum value over the whole imaging period, and 
                     the time at which it occurred, for each individual. 
                     (uses byIndv4Intvl_ValueCalc)}
             }
}
\usage{
traitExtractFeatures(data, individuals = "Snapshot.ID.Tag", times = "DAP", 
                     starts.intvl = NULL, stops.intvl = NULL, sep.intvl = "to", 
                     responses.singletimes = NULL, 
                     responses.rates = NULL, rates.method = "differences", 
                     growth.rates = NULL, suffices.rates = NULL, 
                     water.use = NULL, responses.water = NULL, 
                     responses.total = NULL, suffix.total = NULL, 
                     responses.max = NULL, 
                     times.whole = NULL, 
                     mergedata = NULL, ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns specified by 
            \code{individuals}, \code{times}, the various \code{responses} 
            arguments and the \code{water.use argument}.}
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
 \item{starts.intvl}{A \code{\link{numeric}} giving the times, in terms of values in 
             \code{times}, that are the initial times for a set of intervals for which 
             \code{growth.rates} and \code{water.use} traits are to be obtained. 
             They may also be used to obtain values for single-time traits.}
 \item{stops.intvl}{A \code{\link{numeric}} giving the times, in terms of values in 
             \code{times}, that are the final times for a set of intervals for which 
             \code{growth.rates} and \code{water.use} traits are to be obtained. 
             They may also be used to obtain values for single-time traits.}
 \item{sep.intvl}{A \code{\link{character}} giving the separator to use in combining 
            a \code{starts.intvl} with a \code{stops.intvl} in constructing the suffix 
            to be appended to an interval trait.}
 \item{responses.singletimes}{A \code{\link{character}} specifying the names of the 
             columns containing \code{responses} for which a column of the values is 
             to be formed for each response for each of the unique values in 
             combined \code{starts.intvl} and \code{stops.intvl}. If no interval 
             \code{responses} are specified, then all of the times for 
             \code{responses.singletimes} can be specified using \code{start.intvl}.  
             leaving \code{stops.intvl} set to \code{NULL}.}
 \item{responses.rates}{A \code{\link{character}} specifying the names of the columns 
             containing responses for which growth rates are to be obtained for 
             the intervals specified by \code{starts.intvl} and \code{stops.intvl}. 
             For \code{rates.method} set to \code{differences}, the growth rates will 
             be computed from the column of the \code{response} values whose name is 
             listed in \code{responses.rates}.  For \code{rates.method} set to 
             \code{derivatives}, the growth rates will be computed from a column with 
             the growth rates computed for each time. The name of the column should be 
             a \code{response} listed in \code{responses.rates} to which is appended an 
             element of \code{suffices.rates}.} 
 \item{rates.method}{A \code{\link{character}} specifying the method to use in 
             calculating the growth rates over an interval for \code{response.smoothed}. 
             The two possibilities are \code{"differences"} and \code{"ratesaverages"}. 
             For \code{differences}, the growth rate for an interval is computed 
             by taking differences between the values of a \code{response} for pairs 
             of \code{times}. For \code{ratesaverage}, the growth rate for an interval 
             is computed by  taking weighted averages of growth rates for times within 
             the interval. That is, \code{differences} operates on the original 
             \code{response} and \code{ratesaverage} operates on the growth rates 
             previously calculated from the \code{response}. The \code{ratesaverage} option is 
             most appropriate when the growth rates are calculated using the 
             derivatives of a fitted curve. The option \code{differences} will be 
             more efficient than \code{ratesaverages} when differencing is being used 
             to calculate growth rates.}
 \item{growth.rates}{A \code{\link{character}} giving the growth rates that are 
             to be obtained for the intervals specified by \code{starts.intvl} and 
             \code{stops.intvl}. It should contain one of both of \code{"AGR"} and 
             \code{"RGR"}.} 
 \item{suffices.rates}{A \code{\link{character}} giving the suffices appended to  
             \code{responses.rates} in constructung the column names for the 
             storing the growth rates specified by \code{growth.rates}. If 
             \code{suffices.rates} is \code{NULL}, then \code{"AGR"} and 
             \code{"RGR"} will be used.} 
 \item{water.use}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the water use values that are to be used 
             in computing the water use traits (WU, WUR, WUI) for the intervals 
             specified by \code{starts.intvl} and \code{stops.intvl}.}
 \item{responses.water}{A \code{\link{character}} giving the names of the columns  
             in \code{data} that are to provide the numerator in calculating a 
             \code{WUI} for the intervals specified using \code{starts.intvl} and 
             \code{stops.intvl}. The denominator will be the values in the column 
             in \code{data} whose name is that given by \code{water.use}. 
             See the \code{Value} section for a description of how 
             \code{responses.water} is incorporated into the names constructed for 
             the water use traits.}
 \item{responses.total}{A \code{\link{character}} specifying the names of the 
             columns containing responses for which a column of the values is 
             to be formed by summing the response for each individual over the 
             whole of the imaging period.}
 \item{suffix.total}{A \code{\link{character}} giving the suffix to be appended to  
             an element of \code{responses.total} in constructing the names of the 
             columns in which results of summing each response specified in 
             \code{responses.total}. If \code{NULL}, the suffix will be contructed 
             by joining the minimum and maximum of the the values in 
             \code{starts.intvl} and \code{stops.intvl}, separated by the value of 
             \code{sep.intvl}.} 
 \item{responses.max}{A \code{\link{character}} specifying the names of the 
             columns containing responses for which columns of the values are 
             to be formed that relate to the maximum of the response for each 
             individual over the whole of the imaging period.}
 \item{times.whole}{A \code{\link{numeric}} giving the starts and stop times of imaging. 
             If \code{NULL}, the start time will be the minimum of \code{starts.intvl} 
             and the stop time will be the maximum of \code{stops.intvl}.}
 \item{mergedata}{A \code{\link{data.frame}} containing a column with the name given 
             in \code{individuals} and for which there is only one row for each value 
             given in this column. In general, it will be that the number of rows in 
             \code{mergedata}  is equal to the number of unique values in the column in 
             \code{data} labelled by the value of \code{individuals}, but this is not 
             mandatory. If \code{mergedata} is not \code{NULL}, the values extracted by 
             \code{traitExtractFeatures} will be \code{\link{merge}}d with it.}
 \item{...}{allows passing of arguments to other functions; not used at present.}
}
\value{A \code{\link{data.frame}} that contains an \code{individuals} column and a 
       column for each extracted trait, in addition to any columns in \code{mergedata}. 
       The number of rows in the \code{\link{data.frame}} will equal the number of 
       unique element of the  \code{individuals} column in \code{data}, except when 
       there are extra values in the \code{individuals} column in \code{data}. If the 
       latter applies, then the number of rows will equal the number of unique 
       values in the combined \code{individuals} columns from \code{mergedata} and 
       \code{data}. 
       
       The names of the columns produced by the function are constructed as follows: 
             \enumerate{
               \item{\bold{single times} -- A name for a single-time trait is formed by appending 
                  a full stop to an element of \code{responses.singletimes}, followed by the value of 
                  \code{times} at which the values were observed.}
               \item{\bold{growth rates for a time interval} --  The name for an interval growth 
                  rate is constructed by concatenating the relevant element of \code{responses.rates}, 
                  \code{growth.rates} and a suffix for the time interval, each separated by a full 
                  stop. The interval suffix is formed by joining its \code{starts.intvl} and 
                  \code{stops.intvl} values, separating them by the value of \code{sep.intvl}.}
               \item{\bold{water use traits for a time interval} -- Construction of the names for 
                  the three water traits begins with the value of \code{water.use}. The rate (WUR) 
                  has either \code{R} or \code{.Rate} added to the value of \code{water.use}, the 
                  latter if the value of \code{water.use} contains a full stop. Similarly the index 
                  (WUI) has either \code{I} or \code{.Index} added to it. The WUI also has the 
                  element of \code{responses.water} used in calculating the WUI prefixed to its 
                  name. All three water use traits have a suffix for the interval appended to their 
                  name. This suffix is contructed by joining its \code{starts.intvl} and 
                  \code{stops.intvl}, separated by the value of \code{sep.intvl}.}
               \item{\bold{whole of imaging period} -- The name for whole-of-imaging total is formed 
                  by combining an element of\code{responses.total}  with \code{suffix interval}, 
                  separating them by a full stop.}
               \item{\bold{maximum} -- The name of the column with the maximum values will be the 
                  result of concatenating the \code{responses.max}, \code{"max"} and 
                  \code{suffix.total}, each separated by a full stop. The name of the column with 
                  the value of \code{times} at which the maximum occurred will be the result of 
                  concatenating the \code{responses.max}, \code{"max"} and the value of \code{times},  
                  each separated by a full stop.}
             }

       The \code{\link{data.frame}} is returned invisibly.
}
\author{Chris Brien}
\seealso{\code{\link{getTimesSubset}}, \code{\link{byIndv4Intvl_GRsAvg}}, 
         \code{\link{byIndv4Intvl_GRsDiff}}, \code{\link{byIndv4Intvl_WaterUse}}, \cr
         \code{\link{byIndv_ValueCalc}}.}
\examples{
 #Load dat
 data(tomato.dat)

 #Define DAP constants 
 DAP.endpts   <- c(18,22,27,33,39,43,51)
 nDAP.endpts <- length(DAP.endpts)
 DAP.starts <- DAP.endpts[-nDAP.endpts]
 DAP.stops   <- DAP.endpts[-1]
 DAP.segs <- list(c(DAP.endpts[1]-1, 39), 
                   c(40, DAP.endpts[nDAP.endpts]))
 #Add PSA rates and smooth PSA, also producing sPSA rates
 tom.dat <- byIndv4Times_SplinesGRs(data = tomato.dat, 
                                    response = "PSA", response.smoothed = "sPSA", 
                                    times = "DAP", rates.method = "differences", 
                                    smoothing.method = "log", 
                                    spline.type = "PS", lambda = 1, 
                                    smoothing.segments = DAP.segs)
  
 #Smooth WU
 tom.dat <- byIndv4Times_SplinesGRs(data = tom.dat, 
                                    response = "WU", response.smoothed = "sWU",
                                    rates.method = "none", 
                                    times = "DAP", 
                                    smoothing.method = "direct", 
                                    spline.type = "PS", lambda = 10^(-0.5), 
                                    smoothing.segments = DAP.segs)
 
 #Extract single-valued traits for each individual
 indv.cols <- c("Snapshot.ID.Tag", "Lane", "Position", "Block", "Cart", "AMF", "Zn")
 indv.dat <- subset(tom.dat, subset = DAP == DAP.endpts[1], 
                    select = indv.cols)
 indv.dat <- traitExtractFeatures(data = tom.dat, 
                                  starts.intvl = DAP.starts, stops.intvl = DAP.stops, 
                                   responses.singletimes = "sPSA", 
                                   responses.rates = "sPSA", 
                                   growth.rates = c("AGR", "RGR"), 
                                   water.use = "sWU", responses.water = "sPSA", 
                                   responses.total = "sWU",
                                   responses.max = "sPSA.AGR",
                                   mergedata = indv.dat)
}
\keyword{hplot}
\keyword{manip}