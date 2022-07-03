\name{RicePrepped.dat}
\docType{data}
\alias{RicePrepped.dat}
\title{Prepped data from an experiment to investigate a rice germplasm panel.}
\description{
The data is the full set of Lanes and Positions from an experiment in a Smarthouse at the Plant Accelerator in Adelaide. It is used in the \code{\link{growthPheno-pkg}} as an executable example to illustrate the use of \code{growthPheno}. The experiment and data collection are described in Al-Tamimi et al. (2016) and the data is derived from the \code{\link{data.frame}} in the file \code{00-raw.0254.dat.rda} that is available from Al-Tamimi et al. (2017); halpf od the unprred data is in \code{\link{RiceRaw.dat}}. 

}
\usage{data(RicePrepped.dat)}
\format{A data.frame containing 14784 observations on 32 variables. The names of the columns in the data.frame are: 

\tabular{rlll}{
 Column \tab Name        \tab Class   \tab Description \cr
\cr
  1 \tab Smarthouse      \tab factor \tab the Smarthouse in which a cart occurs. \cr
  2 \tab Snapshot.ID.Tag \tab character \tab a unique identifier for each cart in the experiment. \cr
  3 \tab xDAP            \tab numeric \tab the numbers of days after planting on which the current \cr
    \tab                 \tab         \tab data was observed. \cr
  4 \tab DAT             \tab factor  \tab the numbers of days after treatment on which the current \cr
    \tab                 \tab         \tab data was observed. \cr
  5 \tab xDAT            \tab numeric \tab the numbers of days after treatment on which the current \cr
    \tab                 \tab         \tab data was observed. \cr
  6 \tab cDAT            \tab numeric \tab a centered numeric covariate for DAT. \cr
  7 \tab DAT.diffs       \tab numeric \tab the number of days between this and the previous \cr
    \tab                 \tab         \tab observations (all one for this experiment). \cr
  8 \tab Lane            \tab factor  \tab the Lane in the 24 Lane x 24 Positions grid. \cr
  9 \tab Position        \tab factor  \tab the Position in the 24 Lane x 24 Positions grid. \cr
 10 \tab cPosn           \tab numeric \tab a centered numeric covaariate for Positions. \cr 
 11 \tab cMainPosn       \tab numeric \tab a centered numeric covaariate for Main plots. \cr
 12 \tab Zone            \tab factor  \tab the Zone of 4 Lanes to which the current cart belonged. \cr
 13 \tab cZone           \tab numeric \tab a centered numeric covariate for Zone. \cr
 14 \tab SHZone          \tab factor  \tab the Zone numbered across the two Smarthouses. \cr
 15 \tab ZLane           \tab factor  \tab the number of the Lane within a Zone. \cr
 16 \tab ZMainunit       \tab factor  \tab the number of the Main plot within a Zone. \cr
 17 \tab Subunit         \tab factor  \tab the number of a Cart within a Main plot. \cr
 18 \tab Reps            \tab numeric \tab the replicate number of each Genotype-Salinity combination. \cr
 19 \tab Genotype        \tab factor  \tab the number assigned to the 298 Genotypes  in the experiment. \cr
 20 \tab Salinity        \tab factor  \tab the Salinity treatment (Control, Salt) allocated to a Cart. \cr
 21 \tab PSA             \tab numeric \tab the Projected shoot area (kpixels). \cr
 22 \tab PSA.AGR         \tab numeric \tab the Absolute Growth Rate for the Projected shoot area (kpixels/day). \cr
 23 \tab PSA.RGR         \tab numeric \tab the Relative Growth Rate for the Projected shoot area (per day). \cr
 24 \tab Tr              \tab numeric \tab the amount of water (g) transpired by a plant. \cr
 25 \tab TrR             \tab numeric \tab the rate of water transpireation (g/day) for a plant. \cr
 26 \tab PSA.TUE         \tab numeric \tab the Transpiration Use Efficiency for PSA (kpixels / day) \cr
    \tab                 \tab         \tab for a plant. \cr
 27 \tab sPSA            \tab numeric \tab the smoothed Projected shoot area (kpixels). \cr
 29 \tab sPSA.AGR        \tab numeric \tab the smoothed Absolute Growth Rate for the Projected shoot \cr
    \tab                 \tab         \tab area (kpixels/day). \cr
 29 \tab sPSA.RGR        \tab numeric \tab the smoothed Relative Growth Rate for the Projected shoot \cr 
    \tab                 \tab         \tab area (per day). \cr
 30 \tab sTr             \tab numeric \tab the smoothed amount of water (g) transpired by a plant. \cr
 31 \tab sTrR            \tab numeric \tab the smoothed rate of water transpireation (g/day) for a plant. \cr
 32 \tab sPSA.TUE        \tab numeric \tab the smoothed Transpiration Use Efficiency for PSA (kpixels / day)  \cr
    \tab                 \tab         \tab for a plant. \cr
}}
\references{
Al-Tamimi, N, Brien, C.J., Oakey, H., Berger, B., Saade, S., Ho, Y. S., Schmockel, S. M., Tester, M. and Negrao, S. (2016) New salinity tolerance loci revealed in rice using high-throughput non-invasive phenotyping. \emph{Nature Communications}, \bold{7}, 13342. Retrieved from \doi{10.1038/ncomms13342}.
}
\source{Al-Tamimi N, Brien C, Oakey H, Berger B, Saade S, Ho YS, Schmockel SM, Tester M, Negrao S. (2017) Data from: Salinity tolerance loci revealed in rice using high-throughput non-invasive phenotyping. Retrieved from: \doi{10.5061/dryad.3118j}.}
\keyword{datasets}