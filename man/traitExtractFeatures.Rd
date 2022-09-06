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
               \item{\bold{growth rates for the imaging period overall} -- the average growth rate (AGR 
                     and/or RGR) over the whole imaging period for each individual. 
                     (uses \code{byIndv4Intvl_GRsDiff} or \code{byIndv4Intvl_GRsAvg})}
               \item{\bold{water use traits for the imaging period overall} -- the total water use (WU), 
                     the water use rate (WUR) and the water use index (WUI) for the whole imaging 
                     period for each individual. (uses \code{byIndv4Intvl_WaterUse})}
               \item{\bold{totals for the imaging period overall} -- the total over the whole imaging 
                     period of a trait for each individual. (uses \code{byIndv4Intvl_ValueCalc})}
               \item{\bold{maximum for the imaging period overall} -- the maximum value over the whole 
                     imaging period, and the time at which it occurred, for each individual. 
                     (uses \code{byIndv4Intvl_ValueCalc})}
             }
}
\usage{
traitExtractFeatures(data, individuals = "Snapshot.ID.Tag", times = "DAP", 
                     starts.intvl = NULL, stops.intvl = NULL, 
                     suffices.intvl = NULL, 
                     responses4intvl.rates = NULL, 
                     growth.rates = NULL, 
                     growth.rates.method = "differences", 
                     suffices.growth.rates = NULL, 
                     water.use4intvl.traits = NULL, 
                     responses4water = NULL, 
                     water.trait.types = c("WU", "WUR", "WUI"), 
                     suffix.water.rate = "R", suffix.water.index = "I", 
                     responses4singletimes = NULL, times.single = NULL, 
                     responses4overall.rates = NULL, 
                     water.use4overall.water = NULL, 
                     responses4overall.water = NULL, 
                     responses4overall.totals = NULL, 
                     responses4overall.max = NULL, 
                     intvl.overall = NULL, suffix.overall = NULL, 
                     sep.times.intvl = "to", sep.suffix.times = ".", 
                     sep.growth.rates = ".", sep.water.traits = "", 
                     mergedata = NULL, ...)
}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns specified by 
            \code{individuals}, \code{times}, the various \code{responses} 
            arguments and the \code{water.use argument}.}
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
 \item{starts.intvl}{A \code{\link{numeric}} giving the times, in terms of values in 
             \code{times}, that are the initial times for a set of intervals for which 
             \code{growth.rates} and \code{water.use} traits are to be obtained. 
             These times may also be used to obtain values for single-time traits 
             (see \code{responses4singletimes}).}
 \item{stops.intvl}{A \code{\link{numeric}} giving the times, in terms of values in 
             \code{times}, that are the end times for a set of intervals for which 
             \code{growth.rates} and \code{water.use} traits are to be obtained. 
             These times may also be used to obtain values for single-time traits 
             (see \code{responses4singletimes}).}
 \item{suffices.intvl}{A \code{\link{character}} giving the suffices for intervals 
            specified using \code{starts.intvl} and \code{stops.intvl}. If \code{NULL}, 
            the suffices are automatically generated using \code{starts.intvl}, 
            \code{stops.intvl} and \code{sep.times.intvl}.}
 \item{responses4intvl.rates}{A \code{\link{character}} specifying the names of the columns 
             containing responses for which growth rates are to be obtained for 
             the intervals specified by \code{starts.intvl} and \code{stops.intvl}. 
             For \code{growth.rates.method} set to \code{differences}, the growth rates will 
             be computed from the column of the \code{response} values whose name is 
             listed in \code{responses4intvl.rates}.  For \code{growth.rates.method} set to 
             \code{derivatives}, the growth rates will be computed from a column with 
             the growth rates computed for each time. The name of the column should be 
             a \code{response} listed in \code{responses4intvl.rates} to which is appended an 
             element of \code{suffices.growth.rates}.} 
 \item{growth.rates.method}{A \code{\link{character}} specifying the method to use  
             in calculating the growth rates over an interval for the responses 
             specified by \code{responses4intvl.rates}. 
             The two possibilities are \code{"differences"} and \code{"ratesaverages"}. 
             For \code{differences}, the growth rate for an interval is computed 
             by taking differences between the values of a response for pairs 
             of \code{times}. For \code{ratesaverage}, the growth rate for an interval 
             is computed by  taking weighted averages of growth rates for times within 
             the interval. That is, \code{differences} operates on the  response and 
             \code{ratesaverage} operates on the growth rates previously calculated 
             from the response, so that the appropriate one of these must be in 
             \code{data}. The \code{ratesaverage} option is most appropriate when the 
             growth rates are calculated using the derivatives of a fitted curve. 
             Note that, for responses for which the AGR has been calculated using 
             differences, both methods will give the same result, but the 
             \code{differences} option will  be more efficient than \code{ratesaverages}.}
 \item{growth.rates}{A \code{\link{character}} specifying which growth rates are 
             to be obtained for the intervals specified by \code{starts.intvl} and 
             \code{stops.intvl}. It should contain one of both of \code{"AGR"} and 
             \code{"RGR"}.} 
 \item{suffices.growth.rates}{A \code{\link{character}} giving the suffices appended to  
             \code{responses4intvl.rates} in constructung the column names for the 
             storing the growth rates specified by \code{growth.rates}. If 
             \code{suffices.growth.rates} is \code{NULL}, then \code{"AGR"} and 
             \code{"RGR"} will be used.} 
 \item{water.use4intvl.traits}{A \code{\link{character}} giving the names of the columns in 
             \code{data} that contain the water use values that are to be used 
             in computing the water use traits (WU, WUR, WUI) for the intervals 
             specified by \code{starts.intvl} and \code{stops.intvl}. If there is 
             only one column name, then the \code{WUI} will be calculated using 
             this name for all column names in \code{responses4water}. If there 
             are several  column names in \code{water.use4intvl.traits}, then there must be 
             either one or the same number of names in \code{responses4water}. 
             If both have same number of names, then the two lists of column names 
             will be processed in parallel, so that a single \code{WUI} will be 
             produced for each pair of \code{responses4water} and \code{water.use4intvl.traits} 
             values.}
 \item{responses4water}{A \code{\link{character}} giving the names of the columns  
             in \code{data} that are to provide the numerator in calculating a 
             \code{WUI} for the intervals specified using \code{starts.intvl} and 
             \code{stops.intvl}. The denominator will be the values in the columns 
             in \code{data} whose names are those given by \code{water.use4intvl.traits}. 
             If there is only one column name in \code{responses4water}, then the 
             \code{WUI} will be calculated using this name for all column names in 
             \code{responses4water}. If there are several column names in 
             \code{responses4water}, then there must be either one or the same 
             number of names in \code{water.use4intvl.traits}. If both have same number of 
             names, then the two lists of column names will be 
             processed in parallel, so that a single \code{WUI} will be produced 
             for each pair of \code{responses4water} and \code{water.use4intvl.traits} values. 
             
             See the \code{Value} section for a description of how 
             \code{responses4water} is incorporated into the names constructed for 
             the water use traits.}
 \item{water.trait.types}{A \code{\link{character}} listing the trait types to compute
             and return. It should be some combination of \code{WU}, \code{WUR} 
             and \code{WUI}. See Details in \code{\link{byIndv4Intvl_WaterUse}} for 
             how each is calculated.}
 \item{suffix.water.rate}{A \code{\link{character}} giving the label to be appended 
             to the value of \code{water.use4intvl.traits} to form the name of the WUR.}
 \item{suffix.water.index}{A \code{\link{character}} giving the label to be appended 
             to the value of \code{water.use4intvl.traits} to form the name of the WUI.}
 \item{responses4singletimes}{A \code{\link{character}} specifying the names of the 
             columns containing \code{responses} for which a column of the values is 
             to be formed for each response for each of the times values specified in 
             \code{times.single}. If  \code{times.single} is \code{NULL}, then the  
             unique values in the combined \code{starts.intvl} and \code{stops.intvl} 
             will be used.}
 \item{times.single}{A \code{\link{numeric}} giving the times of imaging, for each of 
             which, the values of each \code{responses4singletimes} will be stored in 
             a column of the resulting \code{\link{data.frame}}. If \code{NULL}, then 
             the  unique values in the combined \code{starts.intvl} and \code{stops.intvl} 
             will be used.}
 \item{responses4overall.rates}{A \code{\link{character}} specifying the names of the  
             columns containing responses for which growth rates are to be obtained 
             for the whole imaging period i.e. the interval specified by 
             \code{intvl.overall}. The settings of \code{growth.rates.method}, 
             \code{growth.rates}, \code{suffices.growth.rates}, 
             \code{sep.growth.rates}, \code{suffix.overall} and \code{intvl.overall} 
             will be used in producing the growth rates. See \code{responses4intvl.rates} 
             for more information about how these arguments are used.}
 \item{water.use4overall.water}{A \code{\link{logical}} indicating whether the overall 
             water.traits are to be obtained. The settings of \code{water.trait.types}, 
             \code{suffix.water.rate}, \code{suffix.water.index}, \code{sep.water.traits}, 
             \code{suffix.overall} and \code{intvl.overall} will be used in producing the 
             overall water traits. See \code{water.use4intvl.traits} for more information about 
             how these arguments are used.}
 \item{responses4overall.water}{A \code{\link{character}} giving the names of the columns  
             in \code{data} that are to provide the numerator in calculating a 
             \code{WUI} for the interval corresponding to the whole imaging period. 
             See \code{response.water} for further details. See \code{responses4water} 
             for more information about how this argument is processed.} 
 \item{responses4overall.totals}{A \code{\link{character}} specifying the names of the 
             columns containing responses for which a column of the values is 
             to be formed by summing the response for each individual over the 
             whole of the imaging period.}
 \item{responses4overall.max}{A \code{\link{character}} specifying the names of the 
             columns containing responses for which columns of the values are 
             to be formed that relate to the maximum of the response for each 
             individual over the whole of the imaging period.}
 \item{intvl.overall}{A \code{\link{numeric}} giving the starts and stop times of imaging. 
             If \code{NULL}, the start time will be the minimum of \code{starts.intvl} 
             and the stop time will be the maximum of \code{stops.intvl}.}
 \item{suffix.overall}{A \code{\link{character}} giving the suffix to be appended to  
             the names of traits that apply to the whole imagng period. It applies to 
             \code{overall.growth.rates}, \code{water.use4overall.water}, 
             \code{responses4overall.water} and \code{responses4overall.totals}. 
             If \code{NULL}, then nothing will be added.} 
 \item{sep.times.intvl}{A \code{\link{character}} giving the separator to use in 
            combining a \code{starts.intvl} with a \code{stops.intvl} in constructing 
            the suffix to be appended to an interval trait. If set to \code{NULL} and 
            there is only one value for each of \code{starts.intvl} and 
            \code{stops.intvl}, then no suffix will be added; otherwise 
            \code{sep.times.intvl} set to \code{NULL} will result in an error.}
 \item{sep.suffix.times}{A \code{\link{character}} giving the separator to use in 
            appending a suffix for times to a trait. For no separator, set to 
            \code{""}.}
 \item{sep.growth.rates}{A \code{\link{character}} giving the character(s) to be 
             used to separate the \code{suffices.growth.rates} value from the 
             \code{responses4intvl.rates} values in constructing the name for a 
             new rate. It is also used for separating \code{responses4water} 
             values from the \code{suffix.water.index}. For no separator, set to 
             \code{""}.}
 \item{sep.water.traits}{A \code{\link{character}} giving the character(s) to be 
             used to separate the \code{suffix.rate} and \code{suffix.index} values 
             from the \code{response} value in constructing the name for a new 
             rate/index. The default of \code{""} results in no separator.}
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
                  a full stop to an element of \code{responses4singletimes}, followed by the value of 
                  \code{times} at which the values were observed.}
               \item{\bold{growth rates for a time interval} --  The name for an interval growth 
                  rate is constructed by concatenating the relevant element of \code{responses4intvl.rates}, 
                  \code{growth.rates} and a suffix for the time interval, each separated by a full 
                  stop. The interval suffix is formed by joining its \code{starts.intvl} and 
                  \code{stops.intvl} values, separating them by the value of \code{sep.times.intvl}.}
               \item{\bold{growth rates for the whole imaging period} --  The name for an interval growth 
                  rate is constructed by concatenating the relevant element of \code{responses4intvl.rates}, 
                  \code{growth.rates} and \code{suffix.overall}, each separated by a full 
                  stop.}
               \item{\bold{water use traits for a time interval} -- Construction of the names for 
                  the three water traits begins with the value of \code{water.use4intvl.traits}. The rate (WUR) 
                  has either \code{R} or the value of \code{suffix.water.rate} added to the value of 
                  \code{water.use4intvl.traits}. Similarly the index (WUI) has either \code{I} or the value of 
                  \code{suffix.water.index} added to it. The WUI also has the element of 
                  \code{responses4water} used in calculating the WUI prefixed to its name. 
                  All three water use traits have a suffix for the interval appended to their names. 
                  This suffix is contructed by joining its \code{starts.intvl} and \code{stops.intvl}, 
                  separated by the value of \code{sep.times.intvl}.}
               \item{\bold{water use traits for the whole imaging period} -- Construction of the names for 
                  the three water traits begins with the value of \code{water.use4intvl.traits}. The rate (WUR) 
                  has either \code{R} or the value of \code{suffix.water.rate} added to the value of 
                  \code{water.use4intvl.traits}. Similarly the index (WUI) has either \code{I} or the value of 
                  \code{suffix.water.index} added to it. The WUI also has the element of 
                  \code{responses4water} used in calculating the WUI prefixed to its name. 
                  All three water use traits have \code{suffix.overall} appended to their names.}
               \item{\bold{the total for the whole of imaging period} -- The name for whole-of-imaging 
                  total is formed by combining an element of\code{responses4overall.totals}  with \code{suffix.overall}, 
                  separating them by a full stop.}
               \item{\bold{maximum for the whole of imaging period} -- The name of the column with the 
                  maximum values will be the result of concatenating the \code{responses4overall.max}, \code{"max"} 
                  and \code{suffix.overall}, each separated by a full stop. The name of the column with 
                  the value of \code{times} at which the maximum occurred will be the result of 
                  concatenating the \code{responses4overall.max}, \code{"max"} and the value of \code{times},  
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
                                  responses4singletimes = "sPSA", 
                                  responses4intvl.rates = "sPSA", 
                                  growth.rates = c("AGR", "RGR"), 
                                  water.use4intvl.traits = "sWU", 
                                  responses4water = "sPSA", 
                                  responses4overall.totals = "sWU",
                                  responses4overall.max = "sPSA.AGR",
                                  mergedata = indv.dat)
}
\keyword{hplot}
\keyword{manip}