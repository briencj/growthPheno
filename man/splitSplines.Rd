\name{splitSplines}
\alias{splitSplines}
\title{Adds the fits, and optionally growth rates computed from derivatives, after 
      fitting splines to a response for an individual stored in a \code{data.frame} 
      in long format}
\description{Uses \code{\link{fitSpline}} to fit a spline to a subset of the values 
             of \code{response} and stores the fitted values in \code{data}. 
             The subsets are those values with the same levels combinations 
             of the \code{\link{factor}}s listed in \code{individuals}. The degree of smoothing 
             is controlled by the tuning parameters \code{df} and \code{lambda}, 
             related to the penalty, and by \code{npspline.segments}. The \code{smoothing.method} 
             provides for \code{direct} and \code{logarithmic} smoothing. 

      The derivatives of the fitted spline can also be obtained, and the 
      Absolute and Relative Growth Rates ( AGR and RGR) computed using them, provided 
      \code{correctBoundaries} is \code{FALSE}. Otherwise, growth rates can be 
      obtained by difference using \code{\link{splitContGRdiff}}.

      The handling of missing values in the observations is controlled via 
      \code{na.x.action} and \code{na.y.action}. If there are not 
      at least four distinct, nonmissing x-values, a warning is issued and 
      all smoothed values and derivatives are set to \code{NA}. 

      The function \code{\link{probeSmoothing}} can be used to investgate the effect 
      the smoothing parameters (\code{smoothing.method}, \code{df} or  
      \code{lambda}) on the smooth that results.

             \strong{Note: this function is soft deprecated and may be removed in 
             future versions. \cr Use \code{\link{byIndv4Times_SplinesGRs}}.}
}
\usage{
splitSplines(data, response, response.smoothed = NULL, x, 
             individuals = "Snapshot.ID.Tag", INDICES = NULL,
             smoothing.method = "direct", smoothing.segments = NULL, 
             spline.type = "NCSS", df=NULL, lambda = NULL, 
             npspline.segments = NULL, 
             correctBoundaries = FALSE, 
             deriv = NULL, suffices.deriv = NULL, extra.rate = NULL, 
             sep = ".", 
             na.x.action="exclude", na.y.action = "exclude", ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{response.smoothed}{A \code{\link{character}} specifying the name of the column 
             containing the values of the smoothed response variable, corresponding 
             to \code{response}. If \code{response.smoothed} is \code{NULL}, then 
             \code{response.smoothed} is set to the \code{response} to which 
             \code{.smooth} is added.}
 \item{x}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the predictor variable.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of \code{response} 
             that correspond to the \code{response} values for an individual 
             (e.g. plant, pot, cart, plot or unit) that are to be smoothed 
             separately. If the columns corresponding to \code{individuals} 
             are not \code{\link{factor}}(s) then they will be coerced to 
             \code{\link{factor}}(s). The subsets are formed 
             using \code{\link{split}}.}
 \item{INDICES}{A pseudonym for \code{individuals}.}
 \item{smoothing.method}{A \code{\link{character}} giving the smoothing method 
      to use. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response} and then 
      back-transforming by taking the exponentional of the fitted values.}
 \item{smoothing.segments}{A named \code{\link{list}}, each of whose components 
      is a numeric pair specifying the first and last values of an 
      \code{x}-interval whose data is to be subjected as an entity to smoothing 
      using splines. The separate smooths will be combined to form a whole 
      smooth for each individual. If \code{smoothing.segments} is \code{NULL}, 
      the data is not segmented for smoothing.}
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
      decreases as \code{npspline.segments} increases. When the data has been 
      segmented for smoothing (\code{smoothing.segments} is not \code{NULL}), 
      an \code{npspline.segments} value can be supplied for each segment.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted spline 
            values are to have the method of Huang (2001) applied 
            to them to correct for estimation bias at the end-points. Note that 
            \code{spline.type} must be \code{NCSS} and \code{lambda} and \code{deriv} 
            must be \code{NULL} for \code{correctBoundaries} to be set to \code{TRUE}.}
 \item{deriv}{A \code{numeric} specifying one or more orders of derivatives 
      that are required.}
 \item{suffices.deriv}{A \code{\link{character}} giving the characters to be 
                       appended to the names of the derivatives. If \code{NULL} 
                       and the derivative is to be retained then \code{smooth.dv} 
                       is appended.}
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
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a \code{\link{factor}} to delimit 
            levels when the levels of \code{individuals} are combined to identify 
            subsets.}
 \item{na.x.action}{A \code{\link{character}} string that specifies the action to 
            be taken when values of \code{x} are \code{NA}. The possible 
            values are \code{fail}, \code{exclude} or \code{omit}. 
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
       \code{response.smoothed}. If \code{deriv} is not \code{NULL}, 
       columns containing the values of the derivative(s) will be added 
       to \code{data}; the name each of these columns will be the value of 
       \code{response.smoothed} with \code{.dvf} appended, where \code{f} 
       is the order of the derivative, or  the value of \code{response.smoothed} 
       with the corresponding element of \code{suffices.deriv} appended. 
       If \code{RGR} is not \code{NULL}, the RGR is calculated as the ratio 
       of value of the first derivative of the fitted spline and the fitted 
       value for the spline. 
       Any pre-existing smoothed and derivative columns in \code{data} will be 
       replaced. The ordering of the \code{data.frame} for the \code{x} 
       values will be preserved as far as is possible; the main difficulty 
       is with the handling of missing values by the function \code{merge}. 
       Thus, if missing values in \code{x} are retained, they will occur at 
       the bottom of each subset of \code{individuals} and the order will be 
       problematic when there are missing values in \code{y} and 
       \code{na.y.action} is set to \code{omit}.
}
\author{Chris Brien}
\seealso{\code{\link{fitSpline}}, \code{\link{probeSmoothing}}, \code{\link{splitContGRdiff}}, 
         \code{\link{smooth.spline}}, \code{\link{predict.smooth.spline}}, \cr
         \code{\link{split}}}
\references{
Eilers, P.H.C and Marx, B.D. (2021) \emph{Practical smoothing: the joys of P-splines}. Cambridge University Press, Cambridge.
 
Huang, C. (2001) Boundary corrected cubic smoothing splines. \emph{Journal of Statistical Computation and Simulation}, \bold{70}, 107-121.}
\examples{
data(exampleData)
#smoothing with growth rates calculated using derivates
longi.dat <- splitSplines(longi.dat, response="PSA", x="xDAP", 
                          individuals = "Snapshot.ID.Tag", 
                          df = 4, deriv=1, suffices.deriv="AGRdv", 
                          extra.rate = c(RGRdv = "RGR"))
#Use P-splines
longi.dat <- splitSplines(longi.dat, response="PSA", x="xDAP", 
                          individuals = "Snapshot.ID.Tag", 
                          spline.type = "PS", lambda = 0.1, npspline.segments = 10, 
                          deriv=1, suffices.deriv="AGRdv", 
                          extra.rate = c(RGRdv = "RGR"))
#with segmented smoothing
longi.dat <- splitSplines(longi.dat, response="PSA", x="xDAP", 
                          individuals = "Snapshot.ID.Tag", 
                          smoothing.segments = list(c(28,34), c(35,42)), df = 5)
}
\keyword{data}
\keyword{manip}