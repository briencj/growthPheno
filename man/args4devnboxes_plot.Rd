\name{args4devnboxes_plot}
\alias{args4devnboxes_plot}
\title{Creates a \code{list} of the values for the options of profile plots for comparing smooths}
\description{Creates a \code{list} of the values for the options of deviations boxplots for comparing smooths.
             Note that \code{plots.by}, \code{facet.x} and \code{facet.y} jointly 
             define the organization of the plots. The default settings are optimized for \code{\link{traitSmooth}} so that, 
             if you want to change any of these from their default settings when using \code{args4devnboxes_plot} with a function 
             other than \code{\link{traitSmooth}}, then it is recommended that you specify all of them to ensure that the complete set 
             has been correctly specified. Otherwise, the default settings will be those shown here and these may be different to the 
             default settings shown for the function with which you are using \code{args4devnboxes_plot}.}
\usage{
args4devnboxes_plot(plots.by = "Type", 
                    facet.x = c("Method", "Tuning"), facet.y = ".", 
                    collapse.facets.x = TRUE, collapse.facets.y = FALSE, 
                    facet.labeller = NULL, facet.scales = "fixed", 
                    angle.x = 0, 
                    which.plots = "none", 
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
             profiles plots and deviations boxplots. The default of \code{"."} results 
             in no split into columns.}
 \item{facet.y}{A \code{\link{character}} giving the \code{\link{factor}}s 
             to be used to form subsets to be plotted in separate rows of the 
             profiles plots and deviations boxplots. The default of \code{"."} results 
             in no split into rows.}
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
 \item{angle.x}{A \code{\link{numeric}} between 0 and 360 that gives the angle of the 
             x-axis text to the x-axis. It can also be set by supplying, in 
             \code{ggplotFuncs}, a \code{theme} function from \code{ggplot2}.}
\item{which.plots}{A \code{\link{logical}} indicating which plots are to be produced. 
             The options are either \code{none} or \code{absolute.deviations} and/or 
             \code{relative.deviations}. Boxplots of the absolute deviations 
             are specified by \code{absolute.boxplots}, the absolute deviations being 
             the values of a trait minus their smoothed values (observed - smoothed). 
             Boxplots of the relative deviations are specified by \code{relative.boxplots}, 
             the relative deviations being the absolute deviations divided by the 
             smoothed values of the trait.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element.
                These functions are applied in creating the \code{ggplot} 
                object for a boxplots.}
 \item{...}{allows arguments to be passed to other functions; not used at present.}
} 
\value{A named \code{\link{list}}.}
\author{Chris Brien}
\seealso{\code{\link{traitSmooth}}, \code{\link{probeSmooths}}, \code{\link{plotSmoothsComparison}} and \code{\link{args4chosen_plot}}.}
\examples{
args4devnboxes_plot(plots.by = "Type", 
                    facet.x = "Tuning", 
                    facet.y = c("Smarthouse", "Treatment.1"), 
                    which.plots = "absolute")
}
\keyword{data}
\keyword{manip}