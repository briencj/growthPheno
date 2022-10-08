\name{args4chosen_plot}
\alias{args4chosen_plot}
\title{Creates a \code{list} of the values for the options of profile plots for the chosen smooth}
\description{Creates a \code{list} of the values for the options of profile plots (and boxplots facets) for comparing smooths.
             Note that \code{plots.by}, \code{facet.x}, \code{facet.y} and \code{include.raw} jointly 
             define the organization of the plots. The default settings are optimized for \code{\link{traitSmooth}}.}
\usage{
args4chosen_plot(plots.by = NULL, 
                 facet.x = ".", facet.y = ".", 
                 include.raw = "no", 
                 collapse.facets.x = FALSE, collapse.facets.y = FALSE, 
                 facet.labeller = NULL, facet.scales = "fixed", 
                 breaks.spacing.x = -2, angle.x = 0, 
                 colour = "black", colour.column = NULL, 
                 colour.values = NULL, alpha = 0.3, 
                 addMediansWhiskers = TRUE,
                 ggplotFuncs = NULL, 
                 ...)}
\arguments{
 \item{plots.by}{A \code{\link{character}} that gives the names of the set of 
             \code{\link{factor}}s by which the data is to be grouped and a separate 
             plot produced for each group. If \code{NULL}, no groups are formed. If a 
             set of \code{\link{factor}}s, such as  \code{Type}, \code{Tuning} and 
             \code{Method}, that uniquely index the combinations of the 
              smoothing-parameter values is specified, then groups are formed for each 
             combination of the levels of the these \code{\link{factor}}s, and a separate plot is 
             produced for each combination.}
 \item{facet.x}{A \code{\link{character}} giving the names of the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate columns of the 
             profiles plots. The default of \code{"."} results 
             in no split into columns.}
 \item{facet.y}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate rows of the 
             profiles plots. The default of \code{"."} results 
             in no split into rows.}
 \item{include.raw}{A \code{\link{character}} indicating whether plots of the 
             raw (unsmoothed) trait, corresponding to the plots of the smoothed traits, 
             are to be included in profile plots. The options are \code{no}, \code{alone}, 
             \code{facet.x}, or \code{facet.y}. That is, the plots of the raw traits are 
             plotted separately or as part of either \code{facet.x} or 
             \code{facet.y}.}
 \item{collapse.facets.x}{A \code{\link{logical}} to indicate whether all variables 
             specified by \code{facets.x} are to be collapsed to a single variable. 
             Note that the smoothing-parameters factors, if present, are always collapsed.}
 \item{collapse.facets.y}{A \code{\link{logical}} to indicate whether all variables 
             specified by \code{facets.y} are to be collapsed to a single variable. 
             Note that the smoothing-parameters factors, if present, are always collapsed.}
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
 \item{colour}{A \code{\link{character}} specifying a single colour to use in 
              drawing the lines for the profiles. If colouring according to the 
              values of a variable is required then use \code{colour.column}.}
 \item{colour.column}{A \code{\link{character}} giving the name of a column 
              in \code{data} over whose values the colours of the lines are to be 
              varied. The colours can be specified using \code{colour.values}.}
 \item{colour.values}{A \code{\link{character}} vector specifying the values of 
              the colours to use in drawing the lines for the profiles.  
              If this is a named vector, then the values will be matched based 
              on the names. If unnamed, values will be matched in order 
              (usually alphabetical) within the limits of the scale.}
 \item{alpha }{A \code{\link{numeric}} specifying the degrees of 
            transparency to be used in plotting the responses. 
            It is a ratio in which the denominator specifies the number of 
            points (or lines) that must be overplotted to give a solid cover.}
 \item{addMediansWhiskers}{A \code{\link{logical}} indicating whether plots over 
             time of the medians and outer whiskers are to be added to the plot. 
             The outer whiskers are related to the whiskers on a box-and-whisker 
             and are defined as the median plus (and minus) 1.5 times the 
             interquartile range (IQR). Points lying outside the whiskers are 
             considered to be potential outliers.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} 
                object for a profiles plot.}
 \item{...}{allows arguments to be passed to other functions; not used at present.}
} 
\value{A named \code{\link{list}}.}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}}, \code{\link{probeSmooths}}, \code{\link{plotSmoothsComparison}} and \code{\link{args4profile_plot}}.}
\examples{
args4chosen_plot(plots.by = "Type", 
                 facet.x = "Tuning", facet.y = c("Smarthouse", "Treatment.1"), 
                 include.raw = "facet.x", 
                 alpha  = 0.4, 
                 colour.column = "Method", 
                 colour.values = c("orange", "olivedrab"))
}
\keyword{data}
\keyword{manip}