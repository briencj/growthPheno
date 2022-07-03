\name{exampleData}
\docType{data}
\alias{exampleData}
\alias{raw.dat}
\alias{longi.dat}
\alias{cart.dat}
\title{A small data set to use in function examples}
\description{
Imaging  data for 20 of the plants that were imaged over 14 days from an experiment in a Smarthouse in the Plant Accelerator. Producing these files is illustrated in the \code{Rice} vignette and the data is used as a small example in the \code{growthPheno} manual.
}
\usage{data(exampleData)}
\format{Three \code{data.frames}: 
   \enumerate{
       \item{\code{raw.dat} (280 rows by 33 columns) that contains the imaging data for 20 plants by 14 imaging days as produced by the image processing software;} 
       \item{\code{longi.dat} (280 rows by 37 columns) that contains a modified version of the imaging data for the 20 plants by 14 imaging days in \code{raw.dat};} 
       \item{\code{cart.dat} (20 rows by 14 columns) that contains data summarizing the growth features of the 20 plants produced from the data in \code{longi.dat}.}
}}
\keyword{datasets}