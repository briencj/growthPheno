\name{plotAnom}
\alias{plotAnom}
\alias{anomPlot}
\title{Identifies anomalous individuals and produces longitudinal plots without them and with just them}
\description{Uses \code{\link{intervalValueCalculate}} and the function 
       \code{\link{anom}} to identify anomalous individuals. The user can 
       elect to print the anomalous individuals, a longitudinal profile plot 
       without the anomalous individuals and/or a longitudinal profile plot 
       with only the anomalous individuals. The plots are produced using 
       \code{ggplot}. The plot can be facetted so that a grid of plots is 
       produced. 

       \emph{Warning: \code{anomPlot} will be deprecated in future versions, its synonym \code{plotAnom} 
       being preferred.}
}
\usage{
plotAnom(data, x="xDays+24.16666667", response="Area.smooth.RGR", 
         individuals="Snapshot.ID.Tag", 
         breaks=seq(12, 36, by=2), vertical.line=NULL, 
         groupsFactor=NULL, lower=NULL, upper=NULL, 
         start.time=NULL, end.time=NULL, times.factor = "Days", 
         suffix.interval=NULL, 
         columns.retained=c("Snapshot.ID.Tag", "Smarthouse", "Lane", 
                            "Position", "Treatment.1", "Genotype.ID"),
         whichPrint=c("anomalous","innerPlot","outerPlot"), na.rm=TRUE, ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data to be tested 
             and plotted.}
 \item{x}{A \code{\link{character}} giving the variable to be plotted on 
          the x-axis.}
 \item{response}{A \code{character} specifying the response variable that 
                 is to be tested and plotted on the y-axis.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{factor}(s) that define the subsets of the \code{data} for 
             which each subset corresponds to the \code{response} value for an 
             individual.}
 \item{breaks}{A \code{\link{numeric}} vector giving the breaks to be plotted on 
               the x-axis scale.}
 \item{vertical.line}{A \code{\link{numeric}} giving position on the x-axis 
                      at which a vertical line is to be drawn. If \code{NULL}, 
                      no line is drawn.}
 \item{groupsFactor}{A \code{\link{factor}} giving the name of a 
             \code{factor} that defines groups of individuals between which the 
             test for anomalous individuals can be varied by setting values for 
             one or more of \code{lower}, \code{upper}, \code{start.time} and 
             \code{end.time} to be \code{NULL}, a single value or a set of values whose number 
             equals the number of levels of \code{groupsFactor}. If \code{NULL} or 
             only a dingle value is supplied, 
             the test is the same for all individuals.}
 \item{lower}{A \code{\link{numeric}} such that values in \code{response} below 
              it are  considered to be anomalous. If \code{NULL}, there is no 
              testing for values below the lower bound.}
 \item{upper}{A \code{\link{numeric}} such that values in \code{response} above 
              it are  considered to be anomalous. If \code{NULL}, there is no
              testing for values above the upper bound.}
 \item{start.time}{A \code{\link{numeric}} giving the start of the time interval, 
                   in terms of a level of \code{times.factor}, 
                   during which testing for anomalous values is to occur. 
                   If \code{NULL}, the interval will start with the first 
                   observation.}
 \item{end.time}{A \code{\link{numeric}} giving the end of the time interval, 
                 in terms of a level of \code{times.factor}, 
                 during which testing for anomalous values is to occur.
                 If \code{NULL}, the interval will end with the last 
                 observation.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels should be numeric values stored as 
             characters.}
 \item{suffix.interval}{A \code{\link{character}} giving the suffix to be 
               appended to \code{response} to form the 
               name of the column containing the calculated values. If it is
               \code{NULL} then nothing will be appended.}
 \item{columns.retained}{A \code{\link{character}} giving the names of the 
             columns in \code{data} that are to be retained in the 
             \code{data.frame} of anomalous individuals.}
 \item{whichPrint}{A \code{\link{character}} indicating what is to be printed. 
                   If \code{anomalous} is included, the \code{columns.retained} 
                   are printed for the anomalous individuals.}
 \item{na.rm}{A \code{\link{logical}} indicating whether \code{NA} values should 
              be stripped before the testing proceeds.}
 \item{... }{allows for arguments to be passed to \code{\link{plotLongitudinal}}.}
}
\value{A \code{\link{list}} with three components: 
       \enumerate{
          \item \code{data}, a data frame resulting from the \code{\link{merge}} 
                of \code{data} and the \code{\link{logical}} identifying whether 
                or not an individual is anomalous; 
          \item \code{innerPlot}, an object of class \code{ggplot} storing the 
                longitudinal plot of the individuals that are not anomalous; 
          \item \code{outerPlot}, an object of class \code{ggplot} storing the 
                longitudinal plot of only the individuals that are anomalous.} 
       The name of the column indicating anomalous individuals will be result of 
       concatenating the \code{response}, \code{\link{anom}} and, if it is not 
       \code{NULL},  \code{suffix.interval}, each separated by a full stop. 
       The \code{ggplot} objects can be plotted using \code{print} and can be 
       modified by adding \code{ggplot} functions before printing. If there are 
       no observations to plot, \code{NULL} will be returned for the plot.}
\author{Chris Brien}
\seealso{\code{\link{anom}}, \code{\link{intervalValueCalculate}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
anomalous <- plotAnom(longi.dat, response="Area.smooth.AGR", 
                      lower=2.5, start.time=40, 
                      x = "xDays+35.42857143", vertical.line=29, 
                      breaks=seq(28, 42, by=2), 
                      whichPrint=c("innerPlot"), 
                      y.title="Area.smooth.AGR")
}
\keyword{hplot}

