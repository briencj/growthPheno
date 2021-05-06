\name{longitudinalPrime}
\alias{longitudinalPrime}
\title{Selects a set variables to be retained in a data frame of longitudinal data}
\description{Forms the prime traits by selecting a subset of the traits in a data.frame of 
  imaging data produced by the Lemna Tec Scanalyzer. The imaging traits to be retained 
  are specified using the \code{traits} and \code{labsCamerasViews} arguments. Some imaging 
  traits are divided by 1000 to convert them from pixels to kilopixels. 
  Also added are factors and explanatory variates that might be of use in an analysis.}
\usage{
longitudinalPrime(data, cartId = "Snapshot.ID.Tag", 
                  imageTimes = "Snapshot.Time.Stamp", 
                  timeAfterStart = "Time.after.Planting..d.",
                  idcolumns = c("Genotype.ID","Treatment.1"), 
                  traits = list(all = c("Area", 
                                        "Boundary.Points.To.Area.Ratio", 
                                        "Caliper.Length", "Compactness", 
                                        "Convex.Hull.Area"), 
                                 side = c("Center.Of.Mass.Y", 
                                          "Max.Dist.Above.Horizon.Line")),
                  labsCamerasViews = list(all = c("SV1", "SV2", "TV"),
                                          side = c("SV1", "SV2")), 
                  smarthouse.lev = NULL, 
                  calcWaterLoss = TRUE, pixelsPERcm)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns specified by \code{cartId}, 
             \code{imageTimes}, \code{timeAfterStart}, \code{idcolumns}, \code{traits} and 
             \code{cameras} along with  the following columns: 
          
            Smarthouse, Lane, Position, Weight.Before, Weight.After, Water.Amount, 
            
            Projected.Shoot.Area..pixels.

            The defaults for the arguments to \code{longitudinalPrime} requires a 
            \code{\link{data.frame}} containing the following columns, although not necessarily 
            in the order given here:

            Smarthouse, Lane, Position, Weight.Before, Weight.After, Water.Amount, 
            
            Projected.Shoot.Area..pixels., Area.SV1, Area.SV2, Area.TV, 
            
            Boundary.Points.To.Area.Ratio.SV1, Boundary.Points.To.Area.Ratio.SV2, 

            Boundary.Points.To.Area.Ratio.TV,  Caliper.Length.SV1, 
            
            Caliper.Length.SV2, Caliper.Length.TV, 

            Compactness.SV1, Compactness.SV2, Compactness.TV, 

            Convex.Hull.Area.SV1, Convex.Hull.Area.SV2, Convex.Hull.Area.TV, 
            
            Center.Of.Mass.Y.SV1, Center.Of.Mass.Y.SV2, 

            Max.Dist.Above.Horizon.Line.SV1, Max.Dist.Above.Horizon.Line.SV2.}
 \item{cartId}{A \code{\link{character}} giving the name of the column that contains
               the unique Id for each cart.}
 \item{imageTimes}{A \code{\link{character}} giving the name of the column that contains
               the time that each cart was imaged.}
 \item{timeAfterStart}{A \code{\link{character}} giving the name of the column that 
               contains the time after some nominated starting time e.g. the number of days 
               after planting.}
 \item{idcolumns}{A \code{\link{character}} vector giving the names of the columns that 
                  identify differences between the plants or carts e.g. Genotype.ID, 
                  Treatment.1, Treatment.2.}
 \item{traits}{A \code{\link{character}} or a \code{\link{list}} whose components are  
                  \code{\link{character}}s. Each \code{\link{character}} gives the names of the 
                  columns for imaging \code{traits} whose values are required for each of the 
                  camera-view combinations given in the corresponding \code{\link{list}} 
                  component of \code{labsCamerasViews}. If \code{labsCamerasViews} or a 
                  component of \code{labsCamerasViews} is \code{NULL}, then the contents of 
                  \code{traits} or the coresponding component of \code{traits} are merely 
                  treated as the names of columns to be retained.}
 \item{labsCamerasViews}{A \code{\link{character}} or a \code{\link{list}} whose components are  
                  \code{\link{character}}s. Each \code{\link{character}} gives the labels of the  
                  camera-view combinations for which is required values of each of the imaging 
                  \code{traits} in the corresponding \code{\link{character}} of \code{traits}. 
                  It is assumed that the camera-view labels are appended to the trait names and 
                  separated from the trait names by a full stop (.). If \code{labsCamerasViews} 
                  or a component of \code{labsCamerasViews} is \code{NULL}, then the contents 
                  of the \code{traits} or the coresponding component of \code{traits} are 
                  merely treated as the names of columns to be retained.}
\item{smarthouse.lev}{A \code{\link{character}} vector giving the levels to use 
            for the Smarthouse factor. If \code{NULL} then the unique values in Smarthouse 
            will be used.}
 \item{calcWaterLoss}{A \code{\link{logical}} indicating whether to calculate the Water.Loss. If it is 
           \code{FALSE}, \code{Water.Before}, \code{Water.After} and \code{Water.Amount} will not be in 
           the returned \code{\link{data.frame}}. They can be copied across by listing them in a component 
           of \code{traits} and set the corresponding component of \code{cameras} to \code{NULL}.}
 \item{pixelsPERcm }{A \code{\link{numeric}} giving the number of pixels per cm 
           for the images. \emph{No longer used.}}
}
\value{A \code{\link{data.frame}} containing the columns specified by \code{cartId}, 
             \code{imageTimes}, \code{timeAfterStart}, \code{idcolumns}, \code{traits} and 
             \code{cameras}. The defaults will result in the following columns:
       
       \enumerate{
         \item Smarthouse: factor with levels for the Smarthouse
         \item Lane: factor for lane number in a smarthouse
         \item Position: factor for east/west position in a lane
         \item Days: factor for the number of Days After Planting (DAP)
         \item \code{cartId}: unique code for each cart 
         \item \code{imageTimes}: time at which an image was taken in POSIXct format
         \item Reps: factor indexing the replicates for each combination of the factors in \code{idcolumns} (calculated)
         \item xPosn: numeric for the Positions within a Lane (calculated)
         \item Hour: hour of the day, to 2 decimal places, at which the image was taken (calculated)
         \item xDays: numeric for the DAP that is centred by subtracting the mean of the unique days (calculated)
         \item \code{idcolumns}: the columns listed in \code{idcolumns} that have been converted to factors
         \item Weight.Before: weight of the pot before watering (only if \code{calcWaterLoss} is \code{TRUE})
         \item Weight.After: weight of the pot after watering (only if \code{calcWaterLoss} is \code{TRUE})
         \item Water.Amount: the weight of the water added (= Water.After - Water.Before)  (calculated)
         \item Water.Loss: the difference between Weight.Before for the current imaging and the Weight.After for the previous imaging (calculated unless \code{calcWaterLoss} is \code{FALSE})
         \item Area: the Projected.Shoot.Area..pixels. divided by 1000 (calculated)
         \item Area.SV1: the Projected.Shoot.Area from Side View 1 divided by 1000 (calculated)
         \item Area.SV2: the Projected.Shoot.Area from Side View 2 divided by 1000 (calculated)
         \item Area.TV: the Projected.Shoot.Area from Top View divided by 1000 (calculated)
         \item Boundary.To.Area.Ratio.SV1 
         \item Boundary.To.Area.Ratio.SV2 
         \item Boundary.To.Area.Ratio.TV 
         \item Caliper.Length.SV1 
         \item Caliper.Length.SV2 
         \item Caliper.Length.TV
         \item Compactness.SV1 from Side View 1
         \item Compactness.SV2 from Side View 2
         \item Compactness.TV: from Top View 
         \item Convex.Hull.Area.SV1: area of Side View 1 Convex Hull divided by 1000 (calculated)
         \item Convex.Hull.Area.SV2: area of Side View 2 Convex Hull divided by 1000 (calculated)
         \item Convex.Hull.TV: Convex.Hull.Area.TV divided by 1000 (calculated)
         \item Center.Of.Mass.Y.SV1: Centre of Mass from Side View 1
         \item Center.Of.Mass.Y.SV2: Centre of Mass from Side View 2
         \item Max.Dist.Above.Horizon.Line.SV1: the Max.Dist.Above.Horizon.Line.SV1 converted to cm using \code{pixelsPERcm} (calculated)
         \item Max.Dist.Above.Horizon.Line.SV2: the Max.Dist.Above.Horizon.Line.SV2 converted to cm using \code{pixelsPERcm} (calculated)
}}
\details{
The columns are copied from data, except for those columns in the list under \bold{Value} that have `(calculated)' appended.
}
\author{Chris Brien}
\examples{
data(exampleData)
longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1)

longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                    traits = list(a = "Area", c = "Compactness"),
                                    labsCamerasViews = list(all = c("SV1", "SV2", "TV"), 
                                                            t = "TV"))

longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                    traits = c("Area.SV1", "Area.SV2", "Area.TV", 
                                               "Compactness.TV"),
                                    labsCamerasViews = NULL)

longiPrime.dat <- longitudinalPrime(data=raw.dat, smarthouse.lev=1, 
                                    calcWaterLoss = FALSE, 
                                    traits = list(img = c("Area", "Compactness"), 
                                                  H20 = c("Weight.Before","Weight.After",
                                                          "Water.Amount")),
                                    labsCamerasViews = list(all = c("SV1", "SV2", "TV"), 
                                                            H2O = NULL))
}
\keyword{data}
\keyword{manip}