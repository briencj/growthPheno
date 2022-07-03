\name{designFactors}
\alias{designFactors}
\title{Adds the factors and covariates for a blocked, split-unit design}
\description{Add the following  \code{\link{factor}}s and covariates to a date frame containing imaging data from the Plant Accelerator: Zone, xZone, SHZone, ZLane, ZMainunit, Subunit and xMainPosn. It checks that the numbers of levels of the \code{\link{factor}}s are consistent with the observed numbers of carts and observations.}
\usage{
designFactors(data, insertName = NULL, designfactorMethod = "LanePosition", 
              nzones = 6, nlanesperzone = 4, 
              nmainunitsperlane = 11, nsubunitspermain = 2)}
\arguments{
 \item{data}{A \code{\link{data.frame}} to which are to be added the design 
        \code{\link{factor}}s and covariates and which must contain the following columns: 

            Smarthouse, Snapshot.ID.Tag, xDAP, and, 

        if \code{designfactorMethod = "LanePosition"}, Lane and Position.}
 \item{insertName}{A \code{\link{character}} giving the name of the column 
       in the \code{data.frame} after which the new \code{\link{factor}}s 
       and covariates are to be inserted. If \code{NULL}, they are added 
       after the last column.} 
 \item{designfactorMethod}{A \code{\link{character}} giving the method to 
       use to obtain the columns for the design \code{\link{factor}}s Zone, 
       ZLane, Mainunit and Subunit. For \code{LanePosition}, it is assumed 
       that (i) Lane can be divided into Zone and ZLane, each with nzones 
       and nlanesperzone levels, respectively, and (ii) Position can be 
       divided into Mainunit and Subunit, each with nmainunitsperlane and 
       nmainunitsperlane levels, respectively. The \code{\link{factor}} 
       SHZone is formed by combining Smarthouse and Zone and ZMainunit 
       is formed by combining ZLane and Mainunit. For \code{StandardOrder}, 
       the \code{\link{factor}}s Zone, ZLane, Mainunit, Subunit are 
       generated in standard order, with the levels of Subunit changing for 
       every observation and the levels of subsequent changing only after all 
       combinations of the levels of the \code{\link{factor}}s to its right 
       have been cycled through.} 
 \item{nzones}{A \code{\link{numeric}} giving the number of zones in a smarthouse.}
 \item{nlanesperzone}{A \code{\link{numeric}} giving the number of lanes in each zone.}
 \item{nmainunitsperlane}{A \code{\link{numeric}} giving the number of mainunits in each lane.}
 \item{nsubunitspermain}{A \code{\link{numeric}} giving the number of subunits in a main plot.}
}
\value{A \code{\link{data.frame}} including the columns:
       
       \enumerate{
         \item Smarthouse: \code{\link{factor}} with levels for the Smarthouse
         \item Zone: \code{\link{factor}} dividing the Lanes into groups, usually of 4 lanes
         \item cZone: numeric corresponding to Zone, centred by subtracting the mean of the unique positions
         \item SHZone: \code{\link{factor}} for the combinations of Smarthouse and Zone
         \item ZLane: \code{\link{factor}} for the lanes within a Zone
         \item ZMainunit: \code{\link{factor}} for the main units within a Zone
         \item Subunit: \code{\link{factor}} for the subunits
         \item cMainPosn: numeric for the main-plot positions within a Lane, centred by subtracting the mean of the unique Positions
         \item cPosn: numeric for the Positions within a Lane, centred by subtracting the mean of the unique Positions}
}
\details{
The \code{\link{factor}}s Zone, ZLane, ZMainunit and Subunit are derived for each Smarthouse based on the values of \code{nzones}, \code{nlanesperzone}, \code{nmainunitsperlane}, \code{nsubunitspermain}, Zone being the blocks in the split-unit design. Thus, the number of carts in each Smarthouse must be the product of these values and the number of observations must be the product of the numbers of smarthouse, carts and imagings for each cart. If this is not the case, it may be able to be achieved by including in \code{data} rows for extra observations that have values for the Snapshot.ID.Tag, Smarthouse, Lane, Position and Time.after.Planting..d. and the remaining columns for these rows have missing values (\code{NA}) Then SHZone is formed by combining Smarthouse and Zone and the covariates cZone, cMainPosn and cPosn calculated. The covariate cZone is calculated from Zone and cMainPosn is formed from the mean of cPosn for each main plot.
}
\author{Chris Brien}
\examples{
data(exampleData)
longi.dat <- prepImageData(data = raw.dat, smarthouse.lev = 1)
longi.dat <- designFactors(data = longi.dat, insertName = "Reps",
                           nzones = 1, nlanesperzone = 1, nmainunitsperlane = 10, 
                           designfactorMethod="StandardOrder")
}
\keyword{data}
\keyword{manip}