\name{args4smoothing}
\alias{args4smoothing}
\title{Creates a \code{list} of the values for the smoothing parameters to be passed to a smoothing function}
\description{Creates a \code{\link{list}} of the values for the smoothing parameters to be 
             passed to a smoothing function. Note that \code{smoothing.methods}, 
             \code{spline.types}, \code{df} and \code{lambdas} are combined to define the 
             set of smooths. The default settings are optimized for \code{\link{traitSmooth}} so that, 
             if you want to change any of these from their default settings when using \code{args4smoothing} with a function 
             other than \code{\link{traitSmooth}}, then it is recommended that you specify all of them to ensure that the complete set 
             has been correctly specified. Otherwise, the default settings will be those shown here and these may be different to the 
             default settings shown for the function with which you are using \code{args4smoothing}.}
\usage{
 args4smoothing(smoothing.methods = "logarithmic", 
                spline.types = c("NCSS","PS"), 
                df = 5:7, 
                lambdas = list(PS = round(10^c(-0.5, 0, 0.5, 1), 
                                          digits = 3)),
                smoothing.segments = NULL, 
                npspline.segments = NULL, 
                na.x.action="exclude", na.y.action = "trimx", 
                external.smooths = NULL, 
                correctBoundaries = FALSE, 
                combinations = "allvalid", 
                ...)}
\arguments{
 \item{smoothing.methods}{A \code{\link{character}} giving the smoothing method 
      to use. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response} and then 
      back-transforming by taking the exponentional of the fitted values.}
 \item{spline.types}{A \code{\link{character}} giving the type of spline  
      to use. Currently, the possibilites are (i) \code{"NCSS"}, for natural 
      cubic smoothing splines, and (ii) \code{"PS"}, for P-splines.}
 \item{df}{A \code{\link{numeric}} with at least one value that specifies, for 
      natural cubic smoothing splines (\code{NCSS}), the desired equivalent 
      numbers of degrees of freedom of the smooths (trace of the smoother matrix). 
      Lower values result in more smoothing. If \code{df = NULL}, the amount of 
      smoothing can be controlled by including a component named \code{NCSS} 
      in the \code{\link{list}} for \code{lambdas}. If \code{df} is \code{NULL} 
      and \code{lambda} does not include a component named \code{NCSS}, 
      then an error is issued.}
 \item{lambdas}{A named \code{\link{list}} or a \code{\link{numeric}} specifying 
      the positive penalties to apply in order to control the amount of smoothing. 
      The amount of smoothing decreases as \code{lamda} decreases. 
      If \code{lambdas} is a \code{\link{list}}, then include a 
      components with \code{lambdas} values and  named for each of the specified 
      values of \code{spline.types} for which \code{lambdas} are to be used. 
      If \code{spline.types} includes \code{PS}, then a component named 
      \code{PS} with at least one numeric value must be present. If a 
      \code{\link{numeric}}, then it will be converted to a \code{\link{list}} 
      with the single component named \code{PS}.}
 \item{smoothing.segments}{A named \code{\link{list}}, each of whose components 
      is a numeric pair specifying the first and last values of an 
      \code{times}-interval whose data is to be subjected as an entity to smoothing 
      using splines. The separate smooths will be combined to form a whole 
      smooth for each individual. If \code{get.rates} includes \code{smoothed} 
      or is \code{TRUE}, \code{rates.method} is \code{differences} and 
      \code{ntimes2span} is 2, the smoothed growth rates will be computed over 
      the set of segments; otherwise, they will be computed within segments. 
      If \code{smoothing.segments} is \code{NULL}, the data is not 
      segmented for smoothing.}
 \item{npspline.segments}{A \code{\link{numeric}} specifying, for P-splines (\code{PS}), 
      the number of equally spaced segments between \code{min(x)} and \code{max(x)}, 
      excluding missing values, to use in constructing the B-spline basis for the 
      spline fitting. If \code{npspline.segments} is NULL, \code{npspline.segments} 
      is set to the maximum of 10 and \code{ceiling((nrow(data)-1)/2)} i.e. there will 
      be at least 10 segments and, for more than 22 \code{times} values, there will be 
      half as many segments as there are \code{times} values. The amount of smoothing 
      decreases as \code{npspline.segments} increases. When the data has been 
      segmented for smoothing (\code{smoothing.segments} is not \code{NULL}), 
      an \code{npspline.segments} value can be supplied for each segment.}
 \item{na.x.action}{A \code{character} string that specifies the action to 
            be taken when values of \code{x} are \code{NA}. The possible 
            values are \code{fail}, \code{exclude} or \code{omit}. 
            For \code{exclude} and \code{omit}, predictions and derivatives 
            will only be obtained for nonmissing values of \code{x}. 
            The difference between these two codes is that for \code{exclude} the returned 
            \code{data.frame} will have as many rows as \code{data}, the 
            missing values have been incorporated. }
 \item{na.y.action}{A \code{character} string that specifies the action to 
            be taken when values of \code{y}, or the \code{response}, are 
            \code{NA}.  The possible values are \code{fail}, \code{exclude}, 
            \code{omit}, \code{allx}, \code{trimx}, \code{ltrimx} or 
            \code{rtrimx}. For  all options, except \code{fail}, missing 
            values in \code{y} will be removed before smoothing. 
            For \code{exclude} and \code{omit}, predictions 
            and derivatives will be obtained only for nonmissing values of 
            \code{x} that do not have missing \code{y} values. Again, the 
            difference between these two is that, only for \code{exclude} 
            will the missing values be incorporated into the 
            returned \code{data.frame}. For \code{allx}, predictions and 
            derivatives will be obtained for all nonmissing \code{x}. 
            For \code{trimx}, they will be obtained for all nonmissing 
            \code{x} between the first and last nonmissing \code{y} values 
            that have been ordered for \code{x}; for \code{ltrimx} and 
            \code{utrimx} either the lower or upper missing \code{y} 
            values, respectively, are trimmed.}
 \item{external.smooths}{A \code{\link{data.frame}} containing the one or more smooths of 
      a \code{response} in the column specified by \code{smoothed.response}. 
      Multiple smoooths should be supplied in long.format with the same columns 
      as the \code{\link{smooths.frame}} \code{data}, except for the 
      smoothing-parameter columns \code{Type}, \code{TunePar}, \code{TuneVal}, 
      \code{Tuning} and \code{Method}. Only those smoothing-parameter columns that 
      are to be used in any of \code{plots.by}, \code{plots.group}, \code{facet.x} 
      and \code{facet.y} should be included with labels appropriate to the 
      \code{external.smooths}. Those smoothing-parameter columns not included in 
      \code{external.smooths} will have columns of \code{"Other"} added to 
      \code{external.smooths}.

      The growth rates will be computed by differencing according to the settings of 
      \code{get.rates} and \code{trait.types} in the function that calls 
      \code{args4smoothing}.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted spline 
             values are to have the method of Huang (2001) applied 
             to them to correct for estimation bias at the end-points. Note that 
             \code{spline.type} must be \code{NCSS} and \code{lambda} and \code{deriv} 
             must be \code{NULL} for \code{correctBoundaries} to be set to \code{TRUE}.}
 \item{combinations}{A \code{\link{character}} specifying how the values of the different 
      smoothing parameters are to be combined to specify the smooths that are 
      to be obtained. The option \code{allvalid} results in a smooth for each of the 
      combinations of the values of \code{smoothing.methods}, \code{spline.types}, 
      \code{df} and \code{lambdas} that are valid; the other smoothing.args 
      will be the same for all smooths.  

      The option \code{parallel} specifies that, if set, each of four 
      smoothing parameters, \code{smoothing.methods}, \code{spline.types}, 
      \code{df} and \code{lambdas}, must have the same number of values and that this 
      number is the number of different smooths to be produced. The values of the 
      parameters in the same position within each parameter collectively specify a 
      single smooth. Because the value of only one of \code{df} and \code{lambdas} 
      must be specified for a smooth, one of these must be set to \code{NA} and the 
      other to the desired value for each smooth. If all values for one of them is 
      \code{NA}, then the argument may be omitted or set to \code{NULL}.

      The option \code{single} is for the specification of a single smooth. This will 
      mean that only one of \code{df} or \code{lambdas} should be set.}
 \item{...}{allows arguments to be passed to other functions; not used at present.}
} 
\value{A named \code{\link{list}}.}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}} and \code{\link{probeSmooths}}.}
\examples{

args4smoothing(smoothing.methods = "direct", 
               spline.types = "NCSS", df = NULL, lambdas = NULL, 
               smoothing.segments = NULL, npspline.segments = NULL, 
               combinations = "allvalid")
args4smoothing(smoothing.methods = c("log","dir","log"), 
               spline.types = c("NCSS","NCSS","PS"), 
               df = c(4,5,NA), lambdas = c(NA,NA,0.36), 
              combinations = "parallel")
args4smoothing(smoothing.methods = "log", 
               spline.types = "PS", df = NULL, 
               lambdas = 0.36, combinations = "single")
}
\keyword{data}
\keyword{manip}