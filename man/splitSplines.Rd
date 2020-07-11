\name{splitSplines}
\alias{splitSplines}
\title{Adds the fits, and optionally growth rates computed from derivatives, after 
      fitting natural cubic smoothing splines to subsets of a response to a \code{data.frame}}
\description{Uses \code{\link{fitSpline}} to fit a spline to a subset of the values 
             of \code{response} and stores the fitted values in \code{data}. 
             The subsets are those values with the same levels combinations 
             of the factors listed in \code{INDICES}. The degree of smoothing 
             is controlled by \code{df} and the \code{smoothing.method} provides 
             for \code{direct} amd \code{logarithmic} smoothing. 

      The derivatives of the fitted spline can also be obtained, and the 
      Absolute and Relative Growth Rates ( AGR and RGR) computed using them, provided 
      \code{correctBoundaries} is \code{FALSE}. Otherwise, growth rates can be 
      obtained by difference using \code{\link{splitContGRdiff}}.

       By default, \code{smooth.spline} will issue an error if there are not 
       at least four distinct x-values. On the other hand, 
       \code{\link{fitSpline}} issues a warning and sets all smoothed values 
       and derivatives to \code{NA}. The handling of missing values in the 
       observations is controlled via \code{na.x.action} and \code{na.y.action}.
}
\usage{
splitSplines(data, response, x, INDICES, df = NULL, smoothing.method = "direct", 
             correctBoundaries = FALSE, 
             deriv = NULL, suffices.deriv=NULL, RGR=NULL, AGR=NULL, sep=".", 
             na.x.action="exclude", na.y.action = "exclude", ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column to be smoothed.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} that is to be smoothed.}
 \item{x}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the predictor variable.}
 \item{INDICES}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of \code{response} 
             that are to be smoothed separately. If the columns corresponding to 
             \code{INDICES} are not \code{\link{factor}}(s) then they will be 
             coerced to \code{\link{factor}}(s). The subsets are formed 
             using \code{\link{split}}.}
 \item{df}{A \code{\link{numeric}} specifying the desired equivalent number of degrees 
      of freedom of the smooth (trace of the smoother matrix). Lower values 
      result in more smoothing. If \code{df = NULL}, ordinary leave-one-out 
      cross-validation is used to determine the amount of smooth.}
 \item{smoothing.method}{A \code{\link{character}} giving the smoothing method 
      to use. The two possibilites are (i) \code{"direct"}, for directly 
      smoothing the observed \code{response}, and (ii) \code{"logarithmic"}, for 
      smoothing the \code{log}-transformed \code{response} and then 
      back-transforming by taking the exponentional of the fitted values.}
 \item{correctBoundaries}{A \code{\link{logical}} indicating whether the fitted 
            spline values are to have the method of Huang (2001) applied 
            to them to correct for estimation bias at the end-points. Note that 
            \code{deriv} must be \code{NULL} for \code{correctBoundaries} to be 
            set to \code{TRUE}.}
 \item{deriv}{A \code{numeric} specifying one or more orders of derivatives 
      that are required.}
 \item{suffices.deriv}{A \code{\link{character}} giving the characters to be 
                       appended to the names of the derivatives. If \code{NULL} 
                       and the derivative is to be retained then \code{smooth.dv} 
                       is appended.}
 \item{RGR}{A \code{\link{character}} giving the character to be appended 
            to the smoothed \code{response} to create the RGR name, 
            but only when \code{smoothing.method} is \code{direct} and the 
            RGR is required. When \code{smoothing.method} is \code{direct} 
            and the RGR is required \code{RGR} must not be \code{NULL} and 
            \code{deriv} must include one so that the first derivative is 
            available for calculating it.
            When \code{smoothing.method} is \code{logarithmic}, the RGR is 
            the backtransformed first derivative and so, to obtain it, merely 
            include \code{1} in \code{deriv} and any suffix for it in 
            \code{suffices.deriv}. Leave \code{RGR} set to \code{NULL}.}
 \item{AGR}{A \code{\link{character}} giving the character to be appended 
            to the smoothed \code{response} to create the AGR name, 
            but only when \code{smoothing.method} is \code{logarithmic} and the 
            AGR is required. When \code{smoothing.method} is \code{logarithmic} 
            and the AGR is required \code{AGR} must not be \code{NULL} and 
            \code{deriv} must include one so that the first derivative is 
            available for calculating it.
            When \code{smoothing.method} is \code{direct}, the AGR is the 
            backtransformed first derivative and so, to obtain it, merely 
            include \code{1} in \code{deriv} and any suffix for it in 
            \code{suffices.deriv}. Leave \code{AGR} set to \code{NULL}.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{INDICES} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a factor to delimit 
            levels when the levels of \code{INDICES} are combined to identify 
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
       \code{response} with \code{.smooth} appended to it. If \code{deriv} is 
       not \code{NULL}, columns containing the values of the derivative(s) 
       will be added to \code{data}; the name each of these columns will 
       be the value of \code{response} with \code{.smooth.dvf} appended, 
       where \code{f} is the order of the derivative, or  the value of 
       \code{response} with \code{.smooth.} and the corresponding 
       element of \code{suffices.deriv} appended. If \code{RGR} is not 
       \code{NULL}, the RGR is calculated as the ratio of value of the first 
       derivative of the fitted spline and the fitted value for the spline. 
       Any pre-existing smoothed and derivative columns in \code{data} will be 
       replaced. The ordering of the \code{data.frame} for the \code{x} 
       values will be preserved as far as is possible; the main difficulty 
       is with the handling of missing values by the function \code{merge}. 
       Thus, if missing values in \code{x} are retained, they will occur at 
       the bottom of each subset of \code{INDICES} and the order will be 
       problematic when there are missing values in \code{y} and 
       \code{na.y.action} is set to \code{omit}.
}
\author{Chris Brien}
\seealso{\code{\link{fitSpline}}, \code{\link{smooth.spline}}, 
         \code{\link{predict.smooth.spline}}, \code{\link{splitContGRdiff}}, \code{\link{split}}}
\references{
Huang, C. (2001). Boundary corrected cubic smoothing splines. \emph{Journal of Statistical Computation and Simulation}, \bold{70}, 107-121.}
\examples{
data(exampleData)
longi.dat <- splitSplines(longi.dat, response="Area", x="xDays", 
                          INDICES = "Snapshot.ID.Tag", 
                          df = 4, deriv=1, suffices.deriv="AGRdv", RGR="RGRdv")
}
\keyword{data}
\keyword{manip}