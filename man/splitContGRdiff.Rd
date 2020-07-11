\name{splitContGRdiff}
\alias{splitContGRdiff}
\title{Adds the growth rates calculated continuously over time for subsets 
of a response to a \code{data.frame}}
\description{Uses \code{\link{AGRdiff}}, \code{\link{PGR}} and 
             \code{\link{RGRdiff}} to calculate growth rates continuously 
             over time for subsets of the values of \code{response} and 
             stores the results in \code{data}. The subsets are those values 
             with the same levels combinations of the factors listed in 
             \code{INDICES}.}
\usage{
splitContGRdiff(data, responses, INDICES,
                which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                times.factor = "Days")}
\arguments{
 \item{data}{A \code{\link{data.frame}} containing the columns for which 
             growth rates are to be calculated.}
 \item{responses}{A \code{\link{character}} giving the names of the columns in 
             \code{data} for which growth rates are to be calculated.}
 \item{INDICES}{A \code{\link{character}} giving the name(s) of the 
             \code{\link{factor}}(s) that define the subsets of \code{response} 
             for which growth rates are to be calculated continuously. 
             If the columns corresponding to 
             \code{INDICES} are not \code{\link{factor}}(s) then they will be 
             coerced to \code{\link{factor}}(s). The subsets are formed 
             using \code{\link{by}}.}
 \item{which.rates}{A \code{\link{character}} giving the growth rates that are 
                    to be calculated. It should be a combination 
                    \code{"AGR"}, \code{"PGR"} and \code{"RGR"}.}
 \item{times.factor}{A \code{\link{character}} giving the name of the column in 
             \code{data} containing the factor for times at which the data was 
             collected. Its levels will be used in calculating growth rates and 
             should be numeric values stored as characters.}
 \item{suffices.rates}{A \code{\link{character}} giving the characters to be 
                    appended to the names of the responses to provide the 
                    names of the columns containing the calculated growth rates.            
                    The order of the suffices in \code{suffices.rates} should 
                    correspond to the order of the elements of 
                    \code{which.rates}. If \code{NULL}, the values of 
                    \code{which.rates} are used.}
}
\value{A \code{\link{data.frame}} containing  \code{data} to which has been 
       added a column for the differences between the \code{times.factor}, 
       if it is not already in \code{data}, and 
       columns with growth rates. The name of the column for \code{times.factor} 
       differences will be the \code{times.factor} with \code{".diff"} appended and, 
       for each of the growth-rate columns will 
       be the value of \code{response} with one of \code{".AGR"}, \code{".PGR"} 
       or \code{"RGR"} or the corresponding value from \code{suffices.GR} 
       appended.}
\author{Chris Brien}
\seealso{\code{\link{fitSpline}}, \code{\link{splitSplines}}}
\examples{
data(exampleData)
longi.dat <- splitContGRdiff(longi.dat, response="Area.smooth", 
                             INDICES = "Snapshot.ID.Tag", which.rates=c("AGR", "RGR"))
}
\keyword{data}
\keyword{manip}