\name{splitValueCalculate}
\alias{splitValueCalculate}
\title{Calculates a single value that is a function of an individual's values for a response}
\description{Splits the values of a response into subsets corresponding 
             individuals and applies a function that calculates a single 
             value to each individual's observations.
             It includes the ability to calculate the observation number that 
             corresponds to the calculated value of the function and the assocated 
             values of a \code{factor} or \code{numeric}.}
\usage{
splitValueCalculate(response, weights=NULL, individuals = "Snapshot.ID.Tag", 
                    FUN = "max", which.obs = FALSE, which.values = NULL, 
                    data, na.rm=TRUE, sep=".", ...)}
\arguments{
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} from which the values of \code{FUN} are to be  
              calculated.}
 \item{weights}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the weights to be supplied as \code{w} 
             to \code{FUN}.}
 \item{individuals}{A \code{\link{character}} giving the name(s) of the 
             \code{factor}(s) that define the subsets of the \code{data} for 
             which each subset corresponds to the \code{response} value for an 
             individual.}
 \item{FUN}{A \code{\link{character}} giving the name of the \code{function} 
             that calculates the value for each subset.}
 \item{which.obs}{A \code{\link{logical}} indicating whether or not to determine 
             the observation number corresponding to the value of the function, 
             in addition to the value of the function itself. For this, \code{FUN} 
             must return an observed value of the \code{reponse}, e.g. \code{max}.}
 \item{which.values}{A \code{\link{character}} giving the name of the 
             \code{\link{factor}} or \code{\link{numeric}} whose levels or values 
             are associated with  the observation number whose \code{response} 
             value corresponds to the value of the function. For this, \code{FUN} 
             must return an observed value of the \code{reponse}, e.g. \code{max}.}
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             function is to be calculated.}
 \item{na.rm}{A \code{\link{logical}}  indicating whether \code{NA} values 
              should be stripped before the calculation proceeds.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a factor to delimit 
            levels when the levels of \code{individuals} are combined to 
            identify subsets.}
 \item{... }{allows for arguments to be passed to \code{FUN}.}
}
\value{A \code{\link{data.frame}}, with the same number of rows as there are 
       \code{individuals}, containing a column for the \code{individuals} and 
       a column with the values of the function for the \code{individuals}. 
       If the \code{FUN} yields observed \code{response} values then it 
       is possible to obtain the origin of the \code{FUN} result, using 
       either or both of \code{which.obs} and \code{which.values}. 
       If \code{which.obs} is \code{TRUE}, a column with the observation 
       numbers corresponding to the values of the function \code{FUN} for the 
       \code{individuals} will be added to the \code{data.frame}. 
       If \code{which.values} is set to the name of a \code{\link{factor}} or 
       a \code{\link{numeric}}, the levels of that \code{\link{factor}} or the 
       values of that \code{\link{numeric}} corresponding to each observation 
       number will be added to the \code{data.frame}.
       
       The name of the column with the values of the function will be formed by  
       concatenating the \code{response} and \code{FUN}, separated by a full stop. 
       If \code{which.obs} is \code{TRUE}, the column name for the obervations 
       numbers will have \code{.obs} added after \code{FUN} into the column name 
       for the function values; if \code{which.values} is specified, 
       the column name for these values will have a full stop followed by 
       \code{which.values} added after \code{FUN} into the column name 
       for the function values.
}
\author{Chris Brien}
\seealso{\code{\link{splitContGRdiff}}, \code{\link{splitSplines}}}
\examples{
data(exampleData)
Area.smooth.max.dat <- splitValueCalculate("Area.smooth", data = longi.dat)
AGR.max.dat <- splitValueCalculate("Area.smooth.AGR", FUN="max", data=longi.dat, 
                                    which.values = "Days", which.obs = TRUE)
}
\keyword{data}
\keyword{manip}