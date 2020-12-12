\name{PVA}
\alias{PVA}
\title{Selects a subset of variables using Principal Variable Analysis (PVA)}
\description{Principal Variable Analysis (PVA) (Cumming and Wooff, 2007) selects a subset 
  from a set of the variables such that the variables in the subset are as 
  uncorrelated as possible, in an effort to ensure that all aspects of the 
  variation in the data are covered. }
\usage{
PVA(obj, ...)}
\arguments{
 \item{obj}{A \code{\link{data.frame}} containing the columns of variables from which the 
             selection is to be made.}
 \item{...}{allows passing of arguments to other functions}
}
\value{A \code{\link{data.frame}} giving the results of the variable selection. 
       It will contain the columns \code{Variable}, \code{Selected}, 
       \code{h.partial}, \code{Added.Propn} and \code{Cumulative.Propn}.
}
\details{\code{PVA} is the generic function for the \code{PVA} method. 
         Use methods("PVA") to get all the methods for the PVA generic.

\code{\link{PVA.data.frame}} is a method for a \code{\link{data.frame}}.

\code{\link{PVA.matrix}} is a method for a \code{\link{matrix}}.
}
\author{Chris Brien}
\seealso{\code{\link{PVA.data.frame}}, \code{\link{PVA.matrix}}, \code{\link{intervalPVA}}, \code{\link{rcontrib}}}
\references{
Cumming, J. A. and D. A. Wooff (2007) Dimension reduction via principal variables. \emph{Computational Statistics
and Data Analysis}, \bold{52}, 550--565.
}
\keyword{data}
\keyword{manip}