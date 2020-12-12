\name{intervalPVA.data.frame}
\alias{intervalPVA.data.frame}
\alias{intervalPVA}
\title{Selects a subset of variables using Principal Variable Analysis (PVA), based on the observed values within a specified time interval}
\description{Principal Variable Analysis (PVA) (Cumming and Wooff, 2007) selects a subset 
  from a set of the variables such that the variables in the subset are as 
  uncorrelated as possible, in an effort to ensure that all aspects of the variation 
  in the data are covered. Here, all observations in a specified time interval are 
  used for calculation the correlations on which the selection is based.}
\usage{
\method{intervalPVA}{data.frame}(obj, responses, times.factor = "Days", start.time, end.time, 
            nvarselect = NULL, p.variance = 1, include = NULL, 
            plot = TRUE, ...)}
\arguments{
 \item{obj}{A \code{\link{data.frame}} containing the columns of variables from which the 
             selection is to be made.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} from which the variables are to be selected.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used to identify the subset and 
             should be numeric values stored as characters.}
 \item{start.time}{A \code{numeric} giving the time,
                   in terms of a level of \code{times.factor}, at which the time interval 
                   begins; observations at this time and up to and including \code{end.time}
                   will be included.}
 \item{end.time}{A \code{numeric} giving the time,
                 in terms of levels of \code{times.factor}, at the end of the interval; 
                 observations after this time will not be included.}
 \item{nvarselect}{A \code{\link{numeric}} specifying the number of variables to be selected, 
                   which includes those listed in \code{include}. If \code{nvarselect = 1}, as 
                   many variables are selected as is need to satisfy \code{p.variance}.}
 \item{p.variance}{A \code{\link{numeric}} specifying the minimum proportion of the variance 
                   that the selected variables must account for,}
 \item{include}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for the variables whose selection is mandatory.}
 \item{plot}{A \code{\link{logical}} indicating whether a plot of the cumulative proportion 
            of the variance explained is to be produced.}
 \item{...}{allows passing of arguments to other functions.}
}
\value{A \code{\link{data.frame}} giving the results of the variable selection. 
       It will contain the columns \code{Variable}, \code{Selected}, 
       \code{h.partial}, \code{Added.Propn} and \code{Cumulative.Propn}.
}
\details{
The variable that is most correlated with the other variables is selected first for inclusion. The partial correlation for each of the remaining variables, given the first selected variable, is calculated and the most correlated of these variables is selects for inclusion next. Then the partial correlations are adjust for the second included variables. This process is repeated until the specified criteria have been satisfied. The possibilities are to:
\enumerate{
  \item the default (\code{nvarselect = NULL} and \code{p.variance = 1}) select all variables in 
        increasing order of amount of information they provide;
  \item select exactly \code{nvarselect} variables;
  \item select just enough variables, up to a maximum of \code{nvarselect} variables, to explain 
        at least \code{p.variance}*100 per cent of the total variance.
}
}
\author{Chris Brien}
\seealso{\code{\link{PVA}}, \code{\link{rcontrib}}}
\references{
Cumming, J. A. and D. A. Wooff (2007) Dimension reduction via principal variables. \emph{Computational Statistics
and Data Analysis}, \bold{52}, 550--565.
}
\examples{
data(exampleData)
responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
               "Density", "Compactness.TV", "Compactness.SV")
results <-  intervalPVA(longi.dat, responses, 
                        start.time = "31", end.time = "31",
                        p.variance=0.9, plot = FALSE)
}
\keyword{data}
\keyword{manip}