\name{byIndv_ValueCalc}
\alias{byIndv_ValueCalc}
\title{Calculates a single value that is a function of an individual's values for a response.}
\description{Applies a \code{function} to calculate a single value from 
             an individual's values for a response in a \code{data.frame} in 
             long format. 
             It includes the ability to calculate the observation number 
             that is closest to the calculated value of the function and 
             the assocated values of a \code{\link{factor}} or \code{numeric}.}
\usage{
byIndv_ValueCalc(data, response, individuals = "Snapshot.ID.Tag", 
                 FUN = "max", which.obs = FALSE, which.values = NULL, 
                 addFUN2name = TRUE, 
                 weights=NULL, na.rm=TRUE, sep=".", ...)}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the column from which the 
             function is to be calculated.}
 \item{response}{A \code{\link{character}} giving the name of the column in 
             \code{data} from which the values of \code{FUN} are to be  
              calculated.}
 \item{individuals}{A \code{\link{character}} giving the name of the 
            \code{\link{factor}} that defines the subsets of the \code{data} 
            for which each subset corresponds to the \code{response} values for 
            an individual (plant/cart/plot/unit).}
 \item{FUN}{A \code{\link{character}} giving the name of the \code{function} 
             that calculates the value for each subset.}
 \item{which.obs}{A \code{\link{logical}} indicating whether or not to determine 
             the observation number corresponding to the observed value that is 
             closest to the value of the function, in addition to the value of 
             the function itself. That is, \code{FUN} need not return an observed 
             value of the \code{reponse}, e.g. \code{quantile}. In the case of 
             multiple observed response values satisfying this condition, 
             the first is returned.}
 \item{which.values}{A \code{\link{character}} giving the name of the 
             \code{\link{factor}} or \code{\link{numeric}} whose values are 
             associated with the \code{response} values and whose value is to be 
             returned for the observation number whose \code{response} value 
             corresponds to the observed value closest to the value of the 
             function. That is, \code{FUN} need not return an observed value of 
             the \code{reponse}, e.g. \code{quantile}. In the case of multiple 
             observed response values satisfying this condition, the value of 
             the \code{which.values} vector for the first of these is 
             returned.}
 \item{addFUN2name}{A \code{\link{logical}} that, if \code{TRUE}, indicates that 
             the \code{FUN} name is to be added to the names of the columns in 
             the \code{\link{data.frame}} returned by \code{byIndv4Intvl_ValueCalc}.}
 \item{weights}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the weights to be supplied as \code{w} 
             to \code{FUN}.}
 \item{na.rm}{A \code{\link{logical}}  indicating whether \code{NA} values 
              should be stripped before the calculation proceeds.}
 \item{sep}{A \code{\link{character}} giving the separator to use when the 
            levels of \code{individuals} are combined. This is needed to avoid 
            using a \code{\link{character}} that occurs in a \code{\link{factor}} to delimit 
            levels when the levels of \code{individuals} are combined to 
            identify subsets.}
 \item{... }{allows for arguments to be passed to \code{FUN}.}
}
\value{A \code{\link{data.frame}}, with the same number of rows as there are 
       \code{individuals}, containing a column for the \code{individuals} and 
       a column with the values of the function for the \code{individuals}. 
       It is also possible to determine observaton numbers or the values of 
       another column in \code{data} for the \code{response} values that are 
       closest to the \code{FUN} results, using either or both of 
       \code{which.obs} and \code{which.values}. If \code{which.obs} is 
       \code{TRUE}, a column with observation numbers is included in the  
       \code{\link{data.frame}}. 
       If \code{which.values} is set to the name of a \code{\link{factor}} or 
       a \code{\link{numeric}},a column containing the levels of that 
       \code{\link{factor}} or the values of that \code{\link{numeric}} is 
       included in the  \code{\link{data.frame}}.
       
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
\seealso{\code{\link{byIndv4Intvl_ValueCalc}}, \code{\link{byIndv4Times_GRsDiff}},
         \code{\link{byIndv4Times_SplinesGRs}}}
\examples{
data(exampleData)
sPSA.max.dat <- byIndv_ValueCalc(data=longi.dat, response = "PSA")
AGR.max.dat <- byIndv_ValueCalc(data=longi.dat, response = "sPSA.AGR", 
                                FUN="max", 
                                which.values = "DAP", which.obs = TRUE)
sPSA.dec1.dat <- byIndv_ValueCalc(data=longi.dat, response = "sPSA", 
                                  FUN="quantile", 
                                  which.values = "DAP", probs = 0.1)
}
\keyword{data}
\keyword{manip}