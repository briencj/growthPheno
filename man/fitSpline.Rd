\name{fitSpline}
\alias{fitSpline}
\title{Fits a spline to a \code{response} in a \code{data.frame}, and growth rates can be 
       computed using derivatives}
\description{Uses \code{smooth.spline} to fit a natural cubic smoothing spline  or \code{JOPS} to fit a 
      P-spline to all the values of \code{response} stored in \code{data}. 

      The amount of smoothing can be controlled by tuning parameters, these being 
      related to the penalty. For a natural cubic smoothing spline, these are 
      \code{df} or \code{lambda} and, for a P-spline, it is \code{lambda}. 
      For a P-spline, \code{npspline.segments} also influences the smoothness of the fit. 
      The \code{smoothing.method} provides for \code{direct} and \code{logarithmic} 
      smoothing. The method of Huang (2001) for correcting the fitted spline for 
      estimation bias at the end-points will be applied when fitting using a 
      natural cubic smoothing spline if \code{correctBoundaries} is \code{TRUE}. 

      The derivatives of the fitted spline can also be obtained, and the 
      Absolute and Relative Growth Rates ( AGR and RGR) computed using them, provided 
      \code{correctBoundaries} is \code{FALSE}. Otherwise, growth rates can be 
      obtained by difference using \code{\link{byIndv4Times_GRsDiff}}.

      The handling of missing values in the observations is controlled via 
      \code{na.x.action} and \code{na.y.action}. If there are not 
      at least four distinct, nonmissing x-values, a warning is issued and 
      all smoothed values and derivatives are set to \code{NA}. 

      The function \code{\link{probeSmoothing}} can be used to investgate the effect 
      the smoothing parameters (\code{smoothing.method} and \code{df} or  
      \code{lambda}) on the smooth that results.
}
\usage{
fitSpline(data, response, response.smoothed, x, 
          smoothing.method = "direct", 
          spline.type = "NCSS", df = NULL, lambda = NULL, 
          npspline.segments = NULL, correctBoundaries = FALSE, 
          deriv = NULL, suffices.deriv = NULL, extra.rate = NULL, 
          na.x.action = "exclude", na.y.action = "trimx", ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response}.}
 \item{x}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the predictor variable.}
 \item{smoothing.method}{A \code{\link{character}} giving the smoothing method 
      to use. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response} and then 
      back-transforming by taking the exponentional of the fitted values.}
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
      the number of equally spaced segments between \code{min(x)} and \code{max(x)}, 
      excluding missing values, to use in constructing the B-spline basis for the 
      spline fitting. If \code{npspline.segments} is NULL, \code{npspline.segments} 
      is set to the maximum of 10 and \code{ceiling((nrow(data)-1)/2)} i.e. there will 
      be at least 10 segments and, for more than 22 \code{x} values, there will be 
      half as many segments as there are \code{x} values. The amount of smoothing 
      decreases as \code{npspline.segments} increases.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted spline 
            values are to have the method of Huang (2001) applied 
            to them to correct for estimation bias at the end-points. Note that 
            \code{spline.type} must be \code{NCSS} and \code{lambda} and \code{deriv} 
            must be \code{NULL} for \code{correctBoundaries} to be set to \code{TRUE}.}
 \item{deriv}{A \code{numeric} specifying one or more orders of derivatives 
      that are required.}
 \item{suffices.deriv}{A \code{\link{character}} giving the characters to be 
            appended to \code{response.method} to construct the names of the derivatives. 
            If \code{NULL} and the derivatives are to be retained, then \code{.dv} 
            followed by the order of the derivative is appended to 
            \code{response.method}}.
 \item{extra.rate}{A named \code{\link{character}} nominating a single growth 
            rate (\code{AGR} or \code{RGR}) to be computed using the first 
            derivative, which one being dependent on the \code{smoothing.method}. 
            The name of this element will used as a suffix to be appended to 
            the response when naming the resulting growth rate (see Examples). 
            If unamed, \code{AGR} or \code{RGR} will be used, as appropriate. 
            Note that, for the \code{smoothing.method} set to \code{direct}, 
            the first derivative is the AGR and so \code{extra.rate} must be set 
            to \code{RGR}, which is computed as the \code{AGR / smoothed response}. 
            For the \code{smoothing.method} set to \code{logarithmic}, 
            the first derivative is the RGR and so \code{extra.rate} must be set 
            to \code{AGR}, which is computed as the \code{RGR * smoothed response}. 
            Make sure that \code{deriv} includes one so that the first derivative 
            is available for calculating the \code{extra.rate}.}
 \item{na.x.action}{A \code{\link{character}} string that specifies the action to 
            be taken when values of \code{x} are \code{NA}. The possible 
            values are \code{fail}, \code{exclude} or \code{omit}. 
            For \code{exclude} and \code{omit}, predictions and derivatives 
            will only be obtained for nonmissing values of \code{x}. 
            The difference between these two codes is that for \code{exclude} 
            the returned \code{data.frame} will have as many rows as 
            \code{data}, the missing values have been incorporated. }
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
\value{A \code{\link{list}} with two components named \code{predictions} and 
       \code{fit.spline}.

       The \code{predictions} component is a \code{data.frame} containing \code{x} 
       and the fitted smooth. The names of the columns will be the value of 
       \code{x} and the value of \code{response.smoothed}. The number of rows in 
       the \code{data.frame} will be equal to the number of pairs that have neither 
       a missing \code{x} or \code{response} and the order of code{x} will be the 
       same as the order in \code{data}. If \code{deriv} is not \code{NULL}, 
       columns containing the values of the derivative(s) will be added to the 
       \code{data.frame}; the name each of these columns will be the value of 
       \code{response.smoothed} with \code{.dvf} appended, where \code{f} is the 
       order of the derivative, or  the value of \code{response.smoothed} and the 
       corresponding element of \code{suffices.deriv} appended. If \code{RGR} is 
       not \code{NULL}, the RGR is calculated as the ratio of value of the first 
       derivative of the fitted spline and the fitted value for the spline.

       The \code{fit.spline} component is a \code{\link{list}} with components
       \describe{
         \item{\code{x}:}{the distinct x values in increasing order;}
         \item{\code{y}:}{the fitted values, with boundary values possibly corrected, and corresponding to \code{x};}
         \item{\code{lev}:}{leverages, the diagonal values of the smoother matrix (NCSS only);}
         \item{\code{lambda}:}{the value of lambda (corresponding to \code{spar} for NCSS - see \code{\link{smooth.spline}});} 
         \item{\code{df}:}{the efective degrees of freedom;}
         \item{\code{npspline.segments}:}{the number of equally spaced segments used for smoothing method set to \code{PS};}
         \item{\code{uncorrected.fit}:}{the object returned by \code{\link{smooth.spline}} for \code{smoothing method} set to \code{NCSS} or by \code{JOPS::psNormal} for \code{PS}.}}
}
\author{Chris Brien}
\seealso{\code{\link{splitSplines}}, \code{\link{probeSmoothing}}, 
         \code{\link{byIndv4Times_GRsDiff}}, \code{\link{smooth.spline}}, \cr
         \code{\link{predict.smooth.spline}}, JOPS.}
\references{
Eilers, P.H.C and Marx, B.D. (2021) \emph{Practical smoothing: the joys of P-splines}. Cambridge University Press, Cambridge.
 
Huang, C. (2001) Boundary corrected cubic smoothing splines. \emph{Journal of Statistical Computation and Simulation}, \bold{70}, 107-121.}
\examples{
data(exampleData)
fit <- fitSpline(longi.dat, response="PSA", response.smoothed = "sPSA", 
                 x="xDAP", df = 4,
                 deriv=c(1,2), suffices.deriv=c("AGRdv","Acc"))
fit <- fitSpline(longi.dat, response="PSA", response.smoothed = "sPSA", 
                 x="xDAP", 
                 spline.type = "PS", lambda = 0.1, npspline.segments = 10, 
                 deriv=c(1,2), suffices.deriv=c("AGRdv","Acc"))
fit <- fitSpline(longi.dat, response="PSA", response.smoothed = "sPSA", 
                 x="xDAP", df = 4,
                 deriv=c(1), suffices.deriv=c("AGRdv"), 
                 extra.rate = c(RGR.dv = "RGR"))
}
\keyword{data}
\keyword{manip}