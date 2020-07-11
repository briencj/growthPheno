\name{calcTimes}
\alias{calcTimes}
\title{Calculates for a set of times, the time intervals after an origin time and the position of each within a time interval}
\description{For the column specified in imageTimes, having converted it to \code{POSIXct} if not already converted, calculates for each value the number of \code{intervalUnits} between the time and the \code{startTime}. Then the number of \code{timePositions} within the \code{intervals} is calculated for the values in \code{imageTimes}. The function \code{difftimes} is used in doing the calculations, but the results are converted to \code{numeric}. For example \code{intervals} could correspond to the number of Days after Planting and the \code{timePositions} to the hour within each day.}
\usage{
calcTimes(data, imageTimes = NULL, timeFormat = "\%Y-\%m-\%d \%H:\%M",
          intervals = "Time.after.Planting..d.", startTime = NULL, 
          intervalUnit = "days", timePositions = NULL)
}
\arguments{
  \item{data}{A \code{\link{data.frame}} containing any columns specified by 
    \code{imageTimes}, \code{intervals} and \code{timePositions}.}
  \item{imageTimes}{A \code{\link{character}} giving the name of the column that contains
    the time that each cart was imaged. Note that in importing data into R, spaces 
    and nonalphanumeric characters in names are converted to full stops. 
    If \code{imageTimes} is \code{NULL} then no calculations are done.}
 \item{timeFormat}{A \code{character} giving the \code{POSIXct} format of characters 
       containing times, in particular \code{imageTimes} and \code{startTime}. 
       Note that if fractions of seconds are required \code{options(digits.secs)} 
       must be used to set the number of decimal places and \code{timeFormat} must 
       use \code{\%OS} for seconds in \code{timeFormat}.} 
  \item{intervals}{A \code{\link{character}} giving the name of the column in \code{data}
    containing, as a \code{\link{numeric}} or a \code{\link{factor}}, 
    the calculated times after \code{startTime} to be 
    plotted on the x-axis.  It is given as the number 
    of \code{intervalUnits} between the two times. If \code{startTime} is 
    \code{NULL} then \code{intervals} is not calculated.}
 \item{startTime}{A \code{character} giving the time, in the \code{POSIXct} 
       format specified by \code{timeFormat}, to be subtracted from 
       \code{imageTimes} to caclualte \code{intervals}. For example, it might be the 
       day of planting or treatment. If \code{startTime} is 
       \code{NULL} then \code{intervals} is not calculated.} 
  \item{intervalUnit}{A \code{\link{character}} giving the name of the unit in which 
    the values of the intervals should be expressed. It must be one of \code{"secs"}, 
    \code{"mins"}, \code{"hours"} or \code{"days"}.}  
  \item{timePositions}{A \code{\link{character}} giving the name of  the column in 
    \code{data} containing, as a \code{\link{numeric}}, the value of the 
    time position within an interval (for example, the time of imaging during the day 
    expressed in hours plus a fraction of an hour). If \code{timePositions} 
    is \code{NULL} then it is not calculated.}
}
\value{A \code{data.frame}, being the unchaged data \code{data.frame} when 
       \code{imageTimes} is \code{NULL} or containing 
       either intervals and/or timePositions depending on which is not \code{NULL}.}
\author{Chris Brien}
\seealso{\code{\link{as.POSIXct}}, \code{\link{imagetimesPlot}}.}
\examples{
  data(exampleData)
  raw.dat <- calcTimes(data = raw.dat, 
                       imageTimes = "Snapshot.Time.Stamp", timePositions = "Hour")
}
\keyword{manip}

