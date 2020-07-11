\name{importExcel}
\alias{importExcel}
\title{Imports an Excel imaging file and allows some renaming of variables}
\description{Uses \code{readxl} to import a sheet of imaging data produced by the 
  Lemna Tec Scanalyzer. Basically, the data consists of imaging data obtained from a 
  set of pots or carts over time. There should be a column, which by default is called 
  \code{Snapshot.ID.Tag}, containing a unique identifier for each cart and a column, 
  which by default is labelled \code{Snapshot.Time.Stamp}, containing 
  the time of imaging for each observation in a row of the sheet. Also, if 
  \code{startTime} is not \code{NULL}, \code{calcTimes} is called to 
  calculate, or recalculate if already present, \code{timeAfterStart} from 
  \code{imageTimes} by subtracting a supplied \code{startTime}. 
  
  Using \code{cameraType}, \code{keepCameraType}, \code{labsCamerasViews} and 
  \code{prefix2suffix}, some flexibility is provided for renaming the columns with 
  imaging data. For example, if the column names are prefixed with 'RGB_SV1', 'RGB_SV2' 
  or 'RGB_TV', the 'RGB_' can be removed and the 'SV1', 'SV2' or 'TV' become suffices. 
  }
\usage{
importExcel(file, sheet="raw data", sep = ",", 
            cartId = "Snapshot.ID.Tag", 
            imageTimes = "Snapshot.Time.Stamp", 
            timeAfterStart = "Time.after.Planting..d.", 
            cameraType = "RGB", keepCameraType = FALSE, 
            labsCamerasViews = NULL, prefix2suffix = TRUE, 
            startTime = NULL,
            timeFormat = "\%Y-\%m-\%d \%H:\%M", 
            plotImagetimes = TRUE, ...)}
\arguments{
 \item{file}{A \code{\link{character}} giving the path and name of the file containing 
             the data.}
 \item{sheet}{A \code{\link{character}} giving the name of the sheet containing 
             the data, that must include columns whose names are as specified by \code{cartId}, 
             which uniquely indexes the carts in the experiment, and 
             \code{imageTimes}, which reflects the time of the imaging from 
             which a particular data value was obtained. It is also assumed that a 
             column whose name is specified by \code{timeAfterStart} is in the sheet or that 
             it will be calculated from \code{imageTimes} using the value of 
             \code{startTime} supplied in the function call.}
 \item{sep}{A \code{\link{character}} giving the separator used in a \code{csv} file.}
 \item{cartId}{A \code{\link{character}} giving the name of the column that contains
               the unique Id for each cart. Note that in importing data into R, spaces 
               and nonalphanumeric characters in names are converted to full stops.}
 \item{imageTimes}{A \code{\link{character}} giving the name of the column that contains
               the time that each cart was imaged. Note that in importing data into R, spaces 
               and nonalphanumeric characters in names are converted to full stops.}
 \item{timeAfterStart}{A \code{\link{character}} giving the name of the column that 
               contains or is to contain the difference between \code{imageTimes} and 
               \code{startTime}. The function \code{calcTimes} is called to calculate the
               differences. For example, it might contain the number of days after 
               planting. Note that in importing data into R, spaces 
               and nonalphanumeric characters in names are converted to full stops.}
 \item{cameraType}{A \code{\link{character}} string nominating the abbreviation used for the 
               cameraType. A warning will be given if no variable names include this cameraType.}
 \item{keepCameraType}{A \code{logical} specifying whether to retain the \code{cameraType} in the 
             variables names. It will be the start of the prefix or suffix and separated from 
             the remander of the prefix or suffix by an underscore (_).}
 \item{labsCamerasViews}{A named \code{\link{character}} whose elements are new labels for the 
               camera-view combinations and the name of each element is the old label for the 
               camera-view combination in the data being imported. If \code{labsCamerasViews} 
               is \code{NULL}, all column names beginning with \code{cameraType} are classed as 
               imaging variables and the unique prefixes amongst them determined. If no imaging 
               variables are found then no changes are made. Note that if you want to include a 
               recognisable \code{cameraType} in a camier-view label, it should be at the start 
               of the the label in \code{labsCamerasViews} and separated from the rest of the 
               label by an underscore (_).}
\item{prefix2suffix}{A \code{logical} specifying whether the variables names with prefixed 
             camera-view labels are to have those prefixes transferred to become suffices. 
             The prefix is assumed to be all the characters up to the first full stop (.) in the 
             variable name and must contain \code{cameraType} to be moved. It is generally 
             assumed that the characters up to the first underscore (_) are the camera type 
             and this is removed if \code{keepCameraType} is \code{FALSE}. If there is no 
             underscore (_), the whole prefix is moved. If \code{labsCamerasViews} is 
             \code{NULL}, all column names beginning with \code{cameraType} are classed as 
             imaging variables and the unique prefixes amongst them determined. 
             If no imaging variables are found then no changes are made.}
 \item{startTime}{A \code{character} giving the time of planting, in the \code{POSIXct} 
       format timeFormat, to be subtracted from \code{imageTimes} 
       in recalculating \code{timeAfterStart}. If \code{startTime} is 
       \code{NULL} then \code{timeAfterStart} is not recalculated.} 
 \item{timeFormat}{A \code{character} giving the \code{POSIXct} format of characters 
       containing times, in particular \code{imageTimes} and \code{startTime}.} 
 \item{plotImagetimes}{A \code{logical} indicating whether a plot of the imaging times against 
             the recalculated \code{Time.After.Planting..d.}. It aids in checking 
             \code{Time.After.Planting..d.} and what occurred in imaging the plants.}
 \item{... }{allows for arguments to be passed to \code{\link{plotImagetimes}}. However, if 
             \code{intervals} is passed an error will occur; use timeAfterStart instead.}
}
\value{A \code{\link{data.frame}} containing the data.}
\author{Chris Brien}
\seealso{\code{\link{as.POSIXct}}, \code{\link{calcTimes}}, \code{\link{plotImagetimes}}}
\examples{\donttest{
filename <- system.file("extdata/rawdata.xlsx", package = "growthPheno", 
                        mustWork = TRUE)
raw.dat <- importExcel(file = filename, 
                       startTime  = "2015-02-11 0:00 AM")
}

camview.labels <- c("SF0", "SL0", "SU0", "TV0")
names(camview.labels) <- c("RGB_Side_Far_0", "RGB_Side_Lower_0", 
                           "RGB_Side_Upper_0", "RGB_TV_0")
filename <- system.file("extdata/raw19datarow.csv", package = "growthPheno", 
                        mustWork = TRUE)
raw.19.dat <- suppressWarnings(importExcel(file = filename,
                                           cartId = "Snapshot.ID.Tags",
                                           startTime = "06/10/2017 0:00 AM",
                                           timeFormat = "\%d/\%m/\%Y \%H:M", 
                                           labsCamerasViews = camview.labels, 
                                           plotImagetimes = FALSE))
}
\keyword{data}
\keyword{manip}