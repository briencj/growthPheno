\name{args4meddevn_plot}
\alias{args4meddevn_plot}
\title{Creates a \code{list} of the values for the options of median deviations plots for smooths}
\description{Creates a \code{list} of the values for the options of median deviations plots for smooths. 
             Note that the arguments \code{plots.by}, \code{plots.group}, \code{facet.x} and \code{facet.y} 
             jointly define the organization of the plots. The default settings are optimized for \code{\link{traitSmooth}} so that, 
             if you want to change any of these from their default settings when using \code{args4meddevn_plot} with a function 
             other than \code{\link{traitSmooth}}, then it is recommended that you specify all of them to ensure that the complete set 
             has been correctly specified. Otherwise, the default settings will be those shown here and these may be different to the 
             default settings shown for the function with which you are using \code{args4meddevn_plot}.
}
\usage{
args4meddevn_plot(plots.by = NULL, plots.group = "Tuning", 
                  facet.x = c("Method","Type"), facet.y = ".",  
                  facet.labeller = NULL, facet.scales = "free_x",
                  breaks.spacing.x = -4, angle.x = 0, 
                  colour.values = NULL, shape.values = NULL, 
                  alpha = 0.5, 
                  propn.note = TRUE, propn.types = NULL, 
                  ggplotFuncs = NULL, 
                  ...)}
\arguments{
 \item{plots.by}{A \code{\link{character}} that give the names of the set of 
             \code{\link{factor}}s by which medians deviations data is to be 
             grouped and a separate plot produced for each group. If \code{NULL}, 
             no groups are formed. If a set of \code{\link{factor}}s, such as  
             \code{Type}, \code{Tuning} and \code{Method}, that uniquely index the 
             combinations of the smoothing-parameter values is specified, then 
             groups are formed for each combination of the levels of the these 
             \code{\link{factor}}s, and a separate plot is produced for each combination.}
 \item{plots.group}{A \code{\link{character}} that gives the names of the set 
             of \code{\link{factor}}s by which the subset of medians deviations data 
             within a single facet in a single plot is to be grouped for plotting  
             as separate lines.}
  \item{facet.x}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate columns of the 
             medians deviations plots. The default of \code{"."} results 
             in no split into columns.}
 \item{facet.y}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate rows of the 
             medians deviations plots. The default of \code{"."} results 
             in no split into rows.}
 \item{facet.labeller}{A \code{\link{ggplot}} \code{\link{function}} for labelling the 
             facets of a plot produced using the \code{\link{ggplot}} function. 
             For more information see \code{\link{ggplot}}.}
 \item{facet.scales}{A \code{\link{character}} specifying whether the scales are shared 
             across all facets of a plot ("fixed"), or do they vary across rows (the 
             default, "free_x"), columns ("free_y"), or both rows and columns ("free")?}
 \item{breaks.spacing.x}{A \code{\link{numeric}} whose absolute values specifies the 
             distance between breaks for the x-axis in a sequence beginning with the 
             minimum x value and continuing up to the maximum x value. If it is negative, 
             the breaks that do not have x values in \code{data} will be omitted. 
             Minor breaks will be at half this value or, if these do not correspond to 
             x-values in \code{data} when \code{breaks.spacing.x} is negative, have 
             a spacing of one. Thus, when \code{breaks.spacing.x} is negative,  
             grid lines will only be included for x-values that occur in \code{data}.  
             These settings can be overwritten by supplying, in \code{ggplotFuncs}, 
             a \code{scale_x_continuous} function from \code{ggplot2}.}
 \item{angle.x}{A \code{\link{numeric}} between 0 and 360 that gives the angle of the 
             x-axis text to the x-axis. It can also be set by supplying, in 
             \code{ggplotFuncs}, a \code{theme} function from \code{ggplot2}.}
 \item{colour.values}{A \code{\link{character}} vector specifying the values of 
              the colours to use in drawing the lines for the medians deviations
              within a facet. If this is a named vector, then the values will be 
              matched based on the names. If unnamed, values will be matched in order 
              (usually alphabetical) within the limits of the scale.}
 \item{shape.values}{A \code{\link{numeric}} vector specifying the values of 
              the shapes to use in drawing the points for the medians deviations
              within a facet. If this is a named vector, then the values will be 
              matched based on the names. If unnamed, values will be matched in 
              order.}
\item{alpha}{A \code{\link{numeric}} specifying the degrees of 
              transparency to be used in plotting a median deviations plot. 
              It is a ratio in which the denominator specifies the number of 
              points (or lines) that must be overplotted to give a solid cover.}
 \item{propn.note}{A \code{\link{logical}} indicating whether a note giving the 
              proportion of the median value of the \code{response} for each time is 
              to be included in the \code{medians.deviations} plots.}
 \item{propn.types}{A \code{\link{numeric}} giving, for each of the 
              \code{trait.types}, the proportion of the median value of the 
              \code{response} for each time to be used to plot envelopes in the 
              median deviations plots. If set to \code{NULL}, the plots of the 
              proprotion envelopes are omitted.}
\item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
              results of evaluating a \code{\link{ggplot}} function. 
              It is created by calling the \code{\link{list}} function with 
              a \code{\link{ggplot}} function call for each element. Note that 
              these functions are applied to the \code{compareian} deviations plots only.}
 \item{...}{allows arguments to be passed to other functions; not used at present.}
} 
\value{A named \code{\link{list}}.}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}}, \code{\link{probeSmooths}} and \code{\link{plotSmoothsMedianDevns}}.}
\examples{
args4meddevn_plot(plots.by = "Type", plots.group = "Tuning",  
                  facet.x = "Method", facet.y = ".",
                  propn.types = c(0.02,0.1, 0.2))
}
\keyword{data}
\keyword{manip}