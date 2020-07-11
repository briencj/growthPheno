\name{twoLevelOpcreate}
\alias{twoLevelOpcreate}
\title{Creates a data.frame formed by applying, for each response, a binary operation 
       to the paired values of two different treatments }
\description{Takes pairs of values for a set of responses indexed by a two-level 
             \code{treatment.factor} and calculates, for each of pair, 
             the result of applying a binary \code{operation} to their values 
             for the two levels of the \code{treatment.factor}. The level of the 
             \code{treatment.factor} designated the \code{control} will be 
             on the right of the binary operator and the value for the other 
             level will be on the left.}
\usage{
twoLevelOpcreate(responses, data, treatment.factor = "Treatment.1", 
                 suffices.treatment = c("Cont","Salt"), control = 1, 
                 columns.suffixed = NULL, 
                 operations = "/", suffices.results="OST", 
                 columns.retained = c("Snapshot.ID.Tag","Smarthouse","Lane", 
                                      "Zones","xZones","SHZones","ZLane",
                                      "ZMainplots","xMainPosn", "Genotype.ID"),
                 by = c("Smarthouse","Zones","ZMainplots"))
}
\arguments{
 \item{responses}{A \code{\link{character}} giving the names of the columns 
             in \code{data} that contain the responses to which the binary 
             operations are to be applied.}
 \item{data}{A \code{\link{data.frame}} containing the columns specified by 
             \code{treatment.factor}, \code{columns.retained} and 
             \code{responses}.}
 \item{treatment.factor}{A \code{\link{factor}} with two levels corresponding to 
                         what is to be designated the \code{control} and treated 
                         observations .}
 \item{suffices.treatment}{A \code{\link{character}} giving the characters to be 
                    appended to the names of the \code{responses} and 
                    \code{columns.suffixed}  in constructing the 
                    names of the columns containing the \code{responses} and 
                    \code{columns.suffixed} for each 
                    level of the \code{treatment.factor}. The order of the 
                    suffices in \code{suffices.treatment} should correspond to 
                    the order of the levels of \code{treatment.factor}.}
 \item{control}{A \code{\link{numeric}}, equal to either 1 or 2, that specifies 
                the level of \code{treatment.factor} that is the control 
                treatment. The value for the \code{control} level will be on 
                the right of the binary operator.}
 \item{columns.suffixed}{A \code{\link{character}} giving the names of the 
             \code{columns.retained} in \code{data} that are to be have the 
             values for each treatment retained and whose names are to be 
             suffixed using \code{suffices.treatment}. Generally, this is done 
             when \code{columns.retained} has different values for different 
             levels of the \code{treatment.factor}.}
 \item{operations}{A \code{\link{character}} giving the binary operations to 
            perform on the values for the two different levels of the 
            \code{treatment.factor}. It should be either of length one, 
            in which case the same operation will be performed for all columns 
            specified in \code{response.GR}, or equal in length to 
            \code{response.GR} so its elements correspond to those of 
            \code{response.GR}.}
 \item{suffices.results}{A \code{\link{character}} giving the characters to be 
                    appended to the names of the responses in constructing the 
                    names of the columns containing the results of applying 
                    the operations. The order of the 
                    suffices in \code{suffices.results} should correspond to 
                    the order of the operators in \code{operations}.}
 \item{columns.retained}{A \code{\link{character}} giving the names of the 
             columns in \code{data} that are to be retained in the 
             \code{data.frame} being created. These are usually 
             \code{\link{factor}s} that index the results of applying the 
             operations and that might be used subsequently.}
 \item{by}{A \code{\link{character}} giving the names of the columns in 
           \code{data} whose combinations will be unique for the observation 
           for each treatment. It is used by \code{\link{merge}} when combining 
           the values of the two treatments in separate 
           columns in the \code{data.frame} to be returned.}
}
\value{A \code{\link{data.frame}} containing the following columns and the values of the :
       
       \enumerate{
         \item those from \code{data} nominated in \code{columns.retained};
         \item those containing the treated values of the columns whose names 
               are specified in \code{responses}; the treated values are 
               those having the other level of \code{treatment.factor} to 
               that specified by \code{control};
         \item those containing the \code{control} values of the columns whose 
               names are specified in \code{responses}; the control values are 
               those having the level of \code{treatment.factor} specified 
               by \code{control};
         \item those containing the values calculated using the binary 
               \code{operations}; the names of these columns will be 
               constructed from \code{responses} by appending 
               \code{suffices.results} to them.
}}
\author{Chris Brien}
\examples{
data(exampleData)
responses <- c("Area.smooth.AGR","Area.smooth.RGR")
cols.retained <-  c("Snapshot.ID.Tag","Smarthouse","Lane","Position",
                    "Days","Snapshot.Time.Stamp", "Hour", "xDays",
                    "Zones","xZones","SHZones","ZLane","ZMainplots",
                    "xMainPosn", "Genotype.ID")
longi.SIIT.dat <- 
  twoLevelOpcreate(responses, longi.dat, suffices.treatment=c("C","S"),
                   operations = c("-", "/"), 
                   suffices.results = c("diff", "SIIT"), 
                   columns.retained = cols.retained, 
                   by = c("Smarthouse","Zones","ZMainplots","Days"))
longi.SIIT.dat <- with(longi.SIIT.dat, 
                            longi.SIIT.dat[order(Smarthouse,Zones,ZMainplots,Days),])
}
\keyword{data}
\keyword{manip}