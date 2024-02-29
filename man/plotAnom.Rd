\name{plotAnom}
\alias{plotAnom}
\title{Identifies anomalous individuals and produces profile plots without them and with just them}
\description{Uses \code{\link{byIndv4Intvl_ValueCalc}} and the function 
       \code{\link{anom}} to identify anomalous individuals in longitudinal 
       data. The user can elect to print the anomalous individuals, a 
       profile plot without the anomalous individuals and/or a profile plot 
       with only the anomalous individuals. The plots are produced using 
       \code{ggplot}. The plot can be facettd so that a grid of plots is 
       produced. 
}
\usage{
plotAnom(data, response="sPSA", 
         individuals="Snapshot.ID.Tag", 
         times = "DAP", x = NULL, 
         breaks.spacing.x = -2, angle.x = 0, 
         vertical.line=NULL, 
         groupsFactor=NULL, lower=NULL, upper=NULL, 
         start.time=NULL, end.time=NULL, 
         suffix.interval=NULL, 
         columns.retained=c("Snapshot.ID.Tag", "Smarthouse", "Lane", 
                            "Position", "Treatment.1", "Genotype.ID"),
         whichPrint=c("anomalous","innerPlot","outerPlot"), na.rm=TRUE, ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the data to be tested 
             and plotted.}
 \item{response}{A \code{character} specifying the response variable that 
                 is to be tested and plotted on the y-axis.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (e.g. plant, pot, cart, plot or unit).}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was collected, 
             either as a \code{\link{numeric}}, \code{\link{factor}}, or 
             \code{\link{character}}. If not a \code{\link{numeric}}, it will be 
             converted to a \code{\link{numeric}} and used to provide the values 
             to be plotted on the x-axis. If a \code{\link{factor}} or 
             \code{\link{character}}, the values should be numerics stored as 
             characters.}
 \item{x}{A \code{\link{character}} specifying a variable, or a function of 
             variables, to be plotted on the x-axis. If \code{NULL}, it will be 
             set to the value of \code{times}, which it can be assumed will be 
             converted to a \code{\link{numeric}}.}
 \item{breaks.spacing.x}{A \code{\link{numeric}} whose absolute values specifies the 
             distance between major breaks for the x-axis in a sequence beginning with the 
             minimum x value and continuing up to the maximum x value. If it is negative, 
             the breaks that do not have x values in \code{data} will be omitted. 
             Minor breaks will be at half major break value or, if these do not correspond to 
             x-values in \code{data} when \code{breaks.spacing.x} is negative, have 
             a spacing of one. Thus, when \code{breaks.spacing.x} is negative,  
             grid lines will only be included for x-values that occur in \code{data}.  
             These settings can be overwritten by supplying, in \code{ggplotFuncs}, 
             a \code{scale_x_continuous} function from \code{ggplot2}.}
 \item{angle.x}{A \code{\link{numeric}} between 0 and 360 that gives the angle of the 
             x-axis text to the x-axis. It can also be set by supplying, in 
             \code{ggplotFuncs}, a \code{theme} function from \code{ggplot2}.}
 \item{vertical.line}{A \code{\link{numeric}} giving position on the x-axis 
                      at which a vertical line is to be drawn. If \code{NULL}, 
                      no line is drawn.}
 \item{groupsFactor}{A \code{\link{factor}} giving the name of a 
             \code{\link{factor}} that defines groups of individuals between which the 
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
                   in terms of a level of \code{times}, 
                   during which testing for anomalous values is to occur. 
                   If \code{NULL}, the interval will start with the first 
                   observation.}
 \item{end.time}{A \code{\link{numeric}} giving the end of the time interval, 
                 in terms of a level of \code{times}, 
                 during which testing for anomalous values is to occur.
                 If \code{NULL}, the interval will end with the last 
                 observation.}
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
                profile plot of the individuals that are not anomalous; 
          \item \code{outerPlot}, an object of class \code{ggplot} storing the 
                profile plot of only the individuals that are anomalous.} 
       The name of the column indicating anomalous individuals will be result of 
       concatenating the \code{response}, \code{\link{anom}} and, if it is not 
       \code{NULL},  \code{suffix.interval}, each separated by a full stop. 
       The \code{ggplot} objects can be plotted using \code{print} and can be 
       modified by adding \code{ggplot} functions before printing. If there are 
       no observations to plot, \code{NULL} will be returned for the plot.}
\author{Chris Brien}
\seealso{\code{\link{anom}}, \code{\link{byIndv4Intvl_ValueCalc}}, \code{\link{ggplot}}.}
\examples{
data(exampleData)
anomalous <- plotAnom(longi.dat, response="sPSA.AGR", 
                      times = "xDAP", 
                      lower=2.5, start.time=40, 
                      vertical.line=29, 
                      breaks.spacing.x = 2, 
                      whichPrint=c("innerPlot"), 
                      y.title="sPSA AGR")
}
\keyword{hplot}

