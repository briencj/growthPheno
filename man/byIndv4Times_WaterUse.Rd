\name{byIndv4Times_WaterUse}
\alias{byIndv4Times_WaterUse}
\title{Adds to a \code{data.frame}, water use traits calculated for consecutive times 
       for individuals in the \code{data.frame} in long format.}
\description{Calculates, for subsets of the data, those water traits out of 
             Water Use (WU), Water Use Rate (WUR) and Water Use Index (WUI) that are 
             specified in \code{which.water.traits}. The traits are calculated 
             continuously over time for the water traits using the values of 
             the variables for consecutive \code{times} stored in \code{data}. 
             The \code{weight.after} column and whichever of  \code{weight.before} 
             and \code{water.added} is not \code{NULL} are used  in the calculations. 
             The subsets are those values with the same levels combinations of the 
             \code{\link{factor}}s listed in \code{individuals}. 
             
             If \code{avail.time.diffs} is \code{FALSE}, the differences between 
             consecutive time values are  calculated for each value of 
            \code{individuals}.}
\usage{
byIndv4Times_WaterUse(data, weight.after = "Weight.After",
                    weight.before = NULL, water.added = NULL, 
                    responseAGR = NULL, 
                    individuals = "Snapshot.ID.Tag", times = "DAP", 
                    which.trait.types = c("WU","WUR","WUI"), 
                    water.trait.names = c("WU", "WUR", "WUI"),  
                    avail.times.diffs = FALSE, ntimes2span = 2)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns to be used in 
             calculating the water traits.}
 \item{weight.after}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the weight of a pot/cart after 
             water has been added.}
 \item{weight.before}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the weight of a pot/cart before 
             water has been added.}
 \item{water.added}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the volume of water that was 
             added to a pot/cart during watering.}
 \item{responseAGR}{A \code{\link{character}} giving the name of the column in 
             \code{data} that contains the values of the Absolute Growth Rate (AGR) 
             of the plant in a a pot/cart in the intervals between successive 
             \code{times} of watering.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
             \code{\link{factor}} that defines the subsets of \code{response} 
             that correspond to the weight, volume and absolute growth rate 
             values for an individual (e.g. plant, pot, cart, plot or unit) for 
             which the water traits are to be calculated continuously. 
             If the column corresponding to \code{individuals} is not a 
             \code{\link{factor}} then it will be coerced to \code{\link{factor}}. 
             The subsets are formed using \code{\link{split}}.}
 \item{times}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the times at which the data was 
             collected, either as a \code{\link{numeric}}, \code{\link{factor}}, 
             or \code{\link{character}}. It will be used in calculating the 
             water traits. If a \code{\link{factor}} or \code{\link{character}}, 
             the values should be numerics stored as characters.}
 \item{which.trait.types}{A \code{\link{character}} giving the water traits that 
             are to be calculated and stored in the \code{\link{data.frame}} 
             \code{data}. It should be a combination of one or more of 
             \code{"WU"}, \code{"WUR"} and \code{"WUI"}.}
 \item{water.trait.names}{A \code{\link{character}} giving the names of the  
             columns in \code{data} for the water traits listed in 
             \code{which.trait.types}.}
 \item{avail.times.diffs}{A \code{\link{logical}} indicating whether there is an 
             appropriate column of \code{times} differences that can be used as 
             the denominator in computing the growth rates. If \code{TRUE}, it will 
             be assumed that the name of the column is the value of  \code{times} 
             with \code{.diffs} appended. If \code{FALSE}, a column, whose 
             column name will be the value of \code{times} with \code{.diffs} 
             appended, will be formed and saved in the result, overwriting any 
             existing columns with the constructed name in \code{data}. It will 
             be calculated using the values of \code{times} in \code{data}.}
 \item{ntimes2span}{A \code{\link{numeric}} giving the number of values in 
             \code{times} to span in calculating growth rates by differencing. 
             Each growth rate is calculated as the difference in the values of 
             one of the \code{responses} for pairs of \code{times} values that 
             are spanned by \code{ntimes2span} \code{times} values divided by 
             the difference between this pair of \code{times} values. For 
             \code{ntimes2span} set to 2, a growth rate is the  
             difference between consecutive pairs of values of one of the 
             \code{responses} divided by the difference between consecutive 
             pairs of \code{times} values.}
}
\value{A \code{\link{data.frame}} containing  \code{data} to which has been 
       added i) a column for the differences between the \code{times}, 
       if it is not already in \code{data}, and (ii) columns with the water 
       traits nominated in \code{which.trait.types} using the names specified in
       \code{water.trait.names}.  The name of the column for \code{times} 
       differences will be the value of \code{times} with \code{".diffs"} 
       appended. Each water trait will be positioned at observation 
       \code{ceiling(ntimes2span + 1) / 2} relative to the two times from 
       which the growth rate is calculated.}

\details{If \code{weight.before} is not \code{NULL}, then the water use (WU)  
         between consecutive \code{times} is calulated as value of 
         \code{weight.after} for first time minus the \code{weight.after} for 
         the second time. If \code{water.added} is not \code{NULL}, then the 
         water use (WU)  between consecutive \code{times} is calulated as value 
         of the \code{water.added} minus the the difference between the 
         \code{weight.after} values for the two  consecutive \code{times}. 
         The \code{WUR} is calculated as the \code{WU} divided by the 
         \code{times.diffs} and the \code{WUI} is calculated as the ratio of
         the values on the column nomiated by \code{responseAGR} to the 
         \code{WUR} values. 
}
\author{Chris Brien}
\seealso{\code{\link{smoothSpline}}, \code{\link{byIndv4Times_GRsDiff}}}
\examples{
data(exampleData)
longi.dat <- byIndv4Times_WaterUse(data = longi.dat, 
                                   weight.before = "Weight.Before", 
                                   individuals = "Snapshot.ID.Tag", 
                                   times = "DAP", 
                                   which.trait.types = c("WU","WUR","WUI"), 
                                   water.trait.names = c("WU", "WUR","WUI"),
                                   responseAGR = "PSA.AGR")
}
\keyword{data}
\keyword{manip}