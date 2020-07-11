\name{designFactors}
\alias{designFactors}
\title{Adds the factors and covariates for a blocked, split-plot design}
\description{Add the following  factors and covariates to a date frame containing imaging data from the Plant Accelerator: Zones, xZones, SHZones, ZLane, ZMainplots, Subplots and xMainPosn. It checks that the numbers of levels of the factors are consistent with the observed numbers of carts and observations.}
\usage{
designFactors(data, insertName = NULL, designfactorMethod = "LanePosition", 
              nzones = 6, nlanesperzone = 4, nmainplotsperlane = 11, nsubplotspermain = 2)}
\arguments{
 \item{data}{A \code{\link{data.frame}} to which are to be added the design factors and covariates and which must contain the following columns: 

    Smarthouse, Snapshot.ID.Tag, XDays, xPosn and, 

    if \code{designfactorMethod = "LanePosition"}, Lane and Position.}
 \item{insertName}{A \code{\link{character}} giving the name of the column in the \code{data.frame} after which the new factors and covariates are to be inserted. If \code{NULL}, they are added after the last column.} 
 \item{designfactorMethod}{A \code{\link{character}} giving the method to use to obtain the columns for the design factors Zones, ZLane, Mainplots and Subplots. For \code{LanePosition}, it is assumed that (i) Lane can be divided into Zones and ZLane, each with nzones and nlanesperzone levels, respectively, and (ii) Position can be divided into Mainplots and Subplots, each with nmainplotsperlane and nmainplotsperlane levels, respectively. The factor SHZones is formed by combining Smarthouse and Zones and ZMainplots is formed by combining ZLane and Mainplots. For \code{StandardOrder}, the factors Zones, ZLane, Mainplots, Subplots are generated in standard order, with the levels of Subplots changing for every observation and the levels of subsequent changing only after all combinations of the levels of the factors to its right have been cycled through.} 
 \item{nzones}{A \code{\link{numeric}} giving the number of zones in a smarthouse.}
 \item{nlanesperzone}{A \code{\link{numeric}} giving the number of lanes in each zone.}
 \item{nmainplotsperlane}{A \code{\link{numeric}} giving the number of mainplots in each lane.}
 \item{nsubplotspermain}{A \code{\link{numeric}} giving the number of subplots in a main plot.}
}
\value{A \code{\link{data.frame}} including the columns:
       
       \enumerate{
         \item Smarthouse: factor with levels for the Smarthouse
         \item Zones: factor dividing the Lanes into groups, usually of 4 lanes
         \item xZones: numeric corresponding to Zones, centred by subtracting the mean of the unique positions
         \item SHZones: factor for the combinations of Smarthouse and Zones
         \item ZLane: factor for the lanes within a Zone
         \item ZMainplots: factor for the main plots within a Zone
         \item Subplots: factor for the subplots
         \item xMainPosn: numeric for the main-plot positions within a Lane, centred by subtracting the mean of the unique positions}
}
\details{
The factors Zones, ZLane, ZMainplots and Subplots are derived for each Smarthouse based on the values of \code{nzones}, \code{nlanesperzone}, \code{nmainplotsperlane}, \code{nsubplotspermain}, Zones being the blocks in the split-plot design. Thus, the number of carts in each Smarthouse must be the product of these values and the number of observations must be the product of the numbers of smarthouse, carts and imagings for each cart. If this is not the case, it may be able to be achieved by including in \code{data} rows for extra observations that have values for the Snapshot.ID.Tag, Smarthouse, Lane, Position and Time.after.Planting..d. and the remaining columns for these rows have missing values (\code{NA}) Then SHZones is formed by combining Smarthouse and Zones and the covariates xZones and xMainPosn calculated. The covariate xZones is calculated from Zones and xMainPosn is formed from the mean of xPosn for each main plot.
}
\author{Chris Brien}
\examples{
data(exampleData)
longi.dat <- designFactors(longi.prime.dat, insertName = "xDays",
                           nzones = 1, nlanesperzone = 1, nmainplotsperlane = 10, 
                           designfactorMethod="StandardOrder")
}
\keyword{data}
\keyword{manip}