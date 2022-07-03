\name{plotImagetimes}
\alias{plotImagetimes}
\title{Plots the position of a time within an interval against the interval for each cart}
\description{Uses \code{ggplot} to produce a plot of the time position within an interval
             against the interval. For example, one might plot the hour of the day carts 
             are imaged against the days after planting (or some other number of 
             days after an event). A line is produced for each value of \code{groupVariable} 
             and the colour is varied according to the value of the \code{colourVariable}. 
             Each \code{Smarthouse} is plotted separately. It aids in checking 
             whether delays  occurred in imaging the plants.
}
\usage{
plotImagetimes(data, intervals = "Time.after.Planting..d.", timePositions = "Hour", 
               groupVariable = "Snapshot.ID.Tag", colourVariable = "Lane", 
               ggplotFuncs = NULL, printPlot = TRUE)
}
\arguments{
  \item{data}{A \code{\link{data.frame}} containing any columns specified by 
    \code{intervals}, \code{timePositions}, \code{groupVariable} and 
    \code{colourVariable}.}
  \item{intervals}{A \code{\link{character}} giving the name of the column in \code{data}
    containing, as a \code{\link{numeric}} or a \code{\link{factor}}, 
    the calculated times to be plotted on the x-axis.  For example, it could be the 
    days after planting or treatment.}
  \item{timePositions}{A \code{\link{character}} giving the name of  the column in 
    \code{data} containing, as a \code{\link{numeric}}, the value of the 
    time position within an interval (for example, the time of imaging during the day 
    expressed in hours plus a fraction of an hour).}
   \item{groupVariable}{A \code{\link{character}} giving the name of  the column in 
                      \code{data} containing the variable to be used to group the 
                      plotting.}
   \item{colourVariable}{A \code{\link{character}} giving the name of  the column 
                      in \code{data} containing the variable to be used to colour 
                      the plotting.}
   \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not to print the 
                  plot.}
}
\value{An object of class "\code{ggplot}", which can be plotted using \code{print}.}
\author{Chris Brien}
\seealso{\code{\link{ggplot}}, \code{\link{calcTimes}}.}
\examples{
data(exampleData)
library(ggplot2)
longi.dat <- calcTimes(longi.dat, imageTimes = "Snapshot.Time.Stamp",
                       timePositions = "Hour")
plotImagetimes(data = longi.dat, intervals = "DAP", timePositions = "Hour",
               ggplotFuncs=list(scale_colour_gradient(low="grey20", high="black"), 
                                geom_line(aes(group=Snapshot.ID.Tag, colour=Lane))))

}
\keyword{hplot}

