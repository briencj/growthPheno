\name{byIndv4Times_SplinesGRs}
\alias{byIndv4Times_SplinesGRs}
\title{For a response in a \code{data.frame} in long format, computes, 
       for a single set of smoothing parameters, smooths of the 
       response, possibly along with growth rates calculated from the 
       smooths.}
\description{Uses \code{\link{smoothSpline}} to fit a spline to the values 
       of \code{response} for each \code{individual} and stores the fitted 
       values in \code{data}. The degree of smoothing is controlled by the 
       tuning parameters \code{df} and \code{lambda}, related to the 
       penalty, and by \code{npspline.segments}. The \code{smoothing.method} 
       provides for \code{direct} and \code{logarithmic} smoothing. 

       The Absolute and Relative Growth Rates ( AGR and RGR) can be computed 
       either using the first derivatives of the splines or by differencing 
       the smooths. If using the first derivative to obtain growth rates, 
       \code{correctBoundaries} must be \code{FALSE}. Derivatives other than the 
       first derivative can also be produced. The function 
       \code{\link{byIndv4Times_GRsDiff}} is used to obtain growth rates by 
       differencing.

       The handling of missing values in the observations is controlled via 
       \code{na.x.action} and \code{na.y.action}. If there are not 
       at least four distinct, nonmissing x-values, a warning is issued and 
       all smoothed values and derivatives are set to \code{NA}. 

       The function \code{\link{probeSmoothing}} can be used to investgate the effect 
       the smoothing parameters \cr (\code{smoothing.method}, \code{df} or  
       \code{lambda}) on the smooth that results.
}
\usage{
byIndv4Times_SplinesGRs(data, response, response.smoothed = NULL, 
                        individuals = "Snapshot.ID.Tag", times, 
                        smoothing.method = "direct", smoothing.segments = NULL, 
                        spline.type = "NCSS", df=NULL, lambda = NULL, 
                        npspline.segments = NULL, 
                        correctBoundaries = FALSE, 
                        rates.method = "differences", 
                        which.rates = c("AGR","RGR"), suffices.rates = NULL, 
                        avail.times.diffs = FALSE, ntimes2span = 2, 
                        extra.derivs = NULL, suffices.extra.derivs=NULL, 
                        sep = ".", 
                        na.x.action="exclude", na.y.action = "trimx", ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response}. If \code{response.smoothed} is \code{NULL}, then 
             \code{response.smoothed} is set to the \code{response} to which is added 
             the prefix \code{s}.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
      \code{\link{factor}}(s) that define the subsets of \code{response} 
      that correspond to the \code{response} values for an individual 
      (plant/cart/plot/unit) that are to be smoothed separately. If the 
      columns corresponding to \code{individuals} are not 
      \code{\link{factor}}(s) then they will be coerced to 
      \code{\link{factor}}(s). The subsets are formed 
      using \code{\link{split}}.}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. It will be used as the values of the predictor 
             variable to be supplied to \code{\link{smooth.spline}} and in 
             calculating growth rates. If a \code{\link{factor}} or 
             \code{\link{character}}, the values should be numerics stored as characters.}
 \item{smoothing.method}{A \code{\link{character}} giving the smoothing method 
      to use. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response} and then 
      back-transforming by taking the exponentional of the fitted values.}
 \item{smoothing.segments}{A named \code{\link{list}}, each of whose components 
      is a numeric pair specifying the first and last values of an 
      \code{times}-interval whose data is to be subjected as an entity to smoothing 
      using splines. The separate smooths will be combined to form a whole 
      smooth for each individual. If \code{get.rates} is \code{TRUE}, 
      \code{rates.method} is \code{differences} and \code{ntimes2span} is 2, 
      the smoothed growth rates will be computed over the set of segments; 
      otherwise, they will be computed within segments. 
      If \code{smoothing.segments} is \code{NULL}, the data is not 
      segmented for smoothing.}
 \item{spline.type}{A \code{\link{character}} giving the type of spline  
      to use. Currently, the possibilites are (i) \code{"NCSS"}, for natural 
      cubic smoothing splines, and (ii) \code{"PS"}, for P-splines.}
 \item{df}{A \code{\link{numeric}} specifying, for natural cubic smoothing splines 
      (\code{NCSS}), the desired equivalent number of degrees of freedom of the 
      smooth (trace of the smoother matrix). Lower values result in more smoothing. 
      If \code{df = NULL}, the amount of smoothing can be controlled by setting 
      \code{lambda}. If both \code{df} and \code{lambda} are \code{NULL}, smoothing 
      is controlled by the default arguments for \code{smooth.spline}, and any 
      that you supply via the ellipsis (\dots) argument.}
 \item{lambda}{A \code{\link{numeric}} specifying the positive penalty to apply. 
      The amount of smoothing decreases as \code{lamda} decreases.}
 \item{npspline.segments}{A \code{\link{numeric}} specifying, for P-splines (\code{PS}), 
      the number of equally spaced segments between \code{min(times)} and \code{max(times)}, 
      excluding missing values, to use in constructing the B-spline basis for the 
      spline fitting. If \code{npspline.segments} is NULL, \code{npspline.segments} 
      is set to the maximum of 10 and \code{ceiling((nrow(data)-1)/2)} i.e. there will 
      be at least 10 segments and, for more than 22 \code{times} values, there will be 
      half as many segments as there are \code{times} values. The amount of smoothing 
      decreases as \code{npspline.segments} increases. When the data has been 
      segmented for smoothing (\code{smoothing.segments} is not \code{NULL}), 
      an \code{npspline.segments} value can be supplied for each segment.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted spline 
            values are to have the method of Huang (2001) applied 
            to them to correct for estimation bias at the end-points. Note that 
            \code{spline.type} must be \code{NCSS} and \code{lambda} and \code{deriv} 
            must be \code{NULL} for \code{correctBoundaries} to be set to \code{TRUE}.}
 \item{rates.method}{A \code{\link{character}} specifying the method to use in 
             calculating the growth rates. The possibilities are 
             \code{none}, \code{differences} and \code{derivatives}.}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
              to be calculated. It should be a combination of one or more of 
             \code{"AGR"}, \code{"PGR"} and \code{"RGR"}.}
 \item{suffices.rates}{A \code{\link{character}} giving the characters to be 
             appended to the names of the responses to provide the names of the 
             columns containing the calculated growth rates. The order of the 
             suffices in \code{suffices.rates} should correspond to the order 
             of the elements of \code{which.rates}. If \code{NULL}, the values 
             of \code{which.rates} are used.}
 \item{avail.times.diffs}{A \code{\link{logical}} indicating whether there is an 
             appropriate column of \code{times} diffserences that can be used as 
             the denominator in computing the growth rates. If \code{TRUE}, it will 
             be assumed that the name of the column is the value of  \code{times} 
             with \code{.diffs} appended. If \code{FALSE}, a column, whose 
             column name will be the value of \code{times} with \code{.diffs} 
             appended, will be formed and saved in the result, overwriting any 
             existing columns with the constructed name in \code{data}. It will 
             be calculated using the values of times in \code{data}.}
 \item{ntimes2span}{A \code{\link{numeric}} giving the number of values in 
             \code{times} to span in calculating growth rates by differencing. 
             Each growth rate is calculated as the difference in the values of 
             one of the \code{responses} for pairs of \code{times} values that 
             are spanned by \code{ntimes2span} \code{times} values divided by 
             the difference between this pair of \code{times} values. For 
             \code{ntimes2span} set to 2, a growth rate is the  
             difference between consecutive pairs of values of one of the 
             \code{responses} divided by the difference between consecutive 
             pairs of \code{times} values.}
 \item{extra.derivs}{A \code{numeric} specifying one or more orders of derivatives 
             that are required, in addition to any required for calculating the growth 
             rates. When \code{rates.method} is \code{derivatives}, these can be 
             derivatives other than the first. Otherwise, any derivatives can be 
             specified.}
 \item{suffices.extra.derivs}{A \code{\link{character}} giving the characters to be 
            appended to \code{response.method} to construct the names of the derivatives. 
            If \code{NULL} and the derivatives are to be retained, then \code{.dv} 
            followed by the order of the derivative is appended to 
            \code{response.method}}.
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a \code{\link{factor}} to delimit 
            levels when the levels of \code{individuals} are combined to identify 
            subsets.}
 \item{na.x.action}{A \code{\link{character}} string that specifies the action to 
            be taken when values of \code{x}, or the \code{times}, are \code{NA}. 
            The possible values are \code{fail}, \code{exclude} or \code{omit}. 
            For \code{exclude} and \code{omit}, predictions and derivatives 
            will only be obtained for nonmissing values of \code{x}. 
            The difference between these two codes is that for \code{exclude} the returned 
            \code{data.frame} will have as many rows as \code{data}, the 
            missing values have been incorporated. }
 \item{na.y.action}{A \code{\link{character}} string that specifies the action to 
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
 \item{...}{allows for arguments to be passed to \code{smooth.spline}.}
}
\value{A \code{\link{data.frame}} containing  \code{data} to which has been 
       added a column with the fitted smooth, the name of the column being 
       the value of \code{response.smoothed}. If \code{rates.method} is 
       not \code{none}, columns for the growth rates listed in 
       \code{which.rates} will be added to \code{data}; the names each of 
       these columns will be the value of \code{response.smoothed} with 
       the elements of \code{which.rates} appended. 
       
       When \code{rates.method} is \code{derivatives} and 
       \code{smoothing.method} is \code{direct}, the \code{AGR} is obtained 
       from the first derivative of the spline for each value of \code{times} 
       and the \code{RGR} is calculated as the \code{AGR} divided by the 
       value of the \code{response.smoothed} for the corresponding time. 
       When \code{rates.method} is \code{derivatives} and 
       \code{smoothing.method} is \code{logarithmic}, the \code{RGR} is 
       obtained from the first derivative of the spline and the \code{AGR} 
       is calculated as the \code{RGR} multiplied by the corresponding 
       value of the \code{response.smoothed}.

       If \code{extra.derivs} is not \code{NULL}, the values for the 
       nominated derivatives will also be added to \code{data}; the names 
       each of these columns will be the value of \code{response.smoothed} 
       with \code{.dvf} appended, where \code{f} is the order of the 
       derivative, or  the value of \code{response.smoothed} 
       with the corresponding element of \code{suffices.deriv} appended. 
       
       Any pre-existing smoothed and growth rate columns in \code{data} will be 
       replaced. The ordering of the \code{data.frame} for the \code{times} 
       values will be preserved as far as is possible; the main difficulty 
       is with the handling of missing values by the function \code{merge}. 
       Thus, if missing values in \code{times} are retained, they will occur at 
       the bottom of each subset of \code{individuals} and the order will be 
       problematic when there are missing values in \code{y} and 
       \code{na.y.action} is set to \code{omit}.
}
\author{Chris Brien}
\seealso{\code{\link{smoothSpline}}, \code{\link{probeSmoothing}}, \code{\link{byIndv4Times_GRsDiff}}, 
         \code{\link{smooth.spline}}, \code{\link{predict.smooth.spline}}, \cr
         \code{\link{split}}}
\references{
Eilers, P.H.C and Marx, B.D. (2021) \emph{Practical smoothing: the joys of P-splines}. Cambridge University Press, Cambridge.
 
Huang, C. (2001) Boundary corrected cubic smoothing splines. \emph{Journal of Statistical Computation and Simulation}, \bold{70}, 107-121.}
\examples{
data(exampleData)
#smoothing with growth rates calculated using derivates
longi.dat <- byIndv4Times_SplinesGRs(data = longi.dat, 
                                     response="PSA", response.smoothed = "sPSA", 
                                     times="DAP", 
                                     df = 4, rates.method = "deriv", 
                                     suffices.rates = c("AGRdv", "RGRdv"))
#Use P-splines
longi.dat <- byIndv4Times_SplinesGRs(data = longi.dat, 
                                     response="PSA", response.smoothed = "sPSA", 
                                     individuals = "Snapshot.ID.Tag", times="DAP", 
                                     spline.type = "PS", lambda = 0.1, 
                                     npspline.segments = 10, 
                                     rates.method = "deriv", 
                                     suffices.rates = c("AGRdv", "RGRdv"))
#with segmented smoothing and no growth rates
longi.dat <- byIndv4Times_SplinesGRs(data = longi.dat, 
                                     response="PSA", response.smoothed = "sPSA", 
                                     individuals = "Snapshot.ID.Tag", times="DAP", 
                                     smoothing.segments = list(c(28,34), c(35,42)), 
                                     df = 5, rates.method = "none")
}
\keyword{data}
\keyword{manip}