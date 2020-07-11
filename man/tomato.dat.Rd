\name{tomato.dat}
\docType{data}
\alias{tomato.dat}
\title{Longitudinal data for an experiment to investigate tomato response to mycorrhizal fungi and zinc}
\description{
The data is from an experiment in a Smarthouse in the Plant Accelerator and is decribed by Watts-Williams et al. (2019). The experiment involves 32 plants, each placed in a pot in a cart, and the carts were assigned 8 treatments using a randomized complete-block design. The main response is Projected Shoot Area (Area for short), being the sum of the plant pixels from three images. The eight treatments were the combinations of 4 Zinc (Zn) levels by two Arbuscular Mycorrhiza Fungi (AMF) levls.  Each plant was imaged on 35 different days after planting (DAPs). It is used to explore the analysis of growth dynamics.

}
\usage{data(tomato.dat)}
\format{A data.frame containing 14784 observations on 33 variables. The names of the columns in the data.frame are: 

\tabular{rlll}{
 Column \tab Name        \tab Class   \tab Description \cr
\cr
  1 \tab Lane            \tab factor \tab the Lane in the 2 Lane x 16 Positions grid. \cr
  2 \tab Position        \tab factor \tab the Position in the 2 Lane x 16 Positions grid. \cr
  3 \tab DAP             \tab factor \tab the numbers of days after planting on which the current \cr
    \tab                 \tab        \tab data was observed. \cr
  4 \tab Snapshot.ID.Tag \tab character \tab a unique identifier for each cart in the experiment. \cr
  5 \tab xDAP            \tab numeric \tab a centered numberic covariate for DAP. \cr
  6 \tab Days.diffs      \tab numeric \tab the number of Days between this and the previous \cr
    \tab                 \tab        \tab observations (all one for this experiment). \cr
  7 \tab xPosn           \tab numeric \tab a centered numeric covaariate for Positions. \cr
  8 \tab Block           \tab factor \tab the block of the randomized complete-block design to \cr
    \tab                 \tab        \tab which the current cart belonged. \cr
  9 \tab Cart            \tab factor \tab the number of the cart within a block. \cr
 10 \tab AMF             \tab factor \tab the AMF treatment (- AMF, +AMF) assigned to the \cr
    \tab                 \tab        \tab cart. \cr
 11 \tab Zn              \tab factor \tab the Zinc level (0, 10, 40, 90) assigned to the cart. \cr
 12 \tab Treatments      \tab factor \tab the combined factor formed from AMF and Zn with\cr 
    \tab                 \tab        \tab levels: (-,0;  -,10; -,40; -,90; +,0;  +,10; +,40; +,90). \cr
 12 \tab Weight.After    \tab numeric \tab the weight of the cart after watering. \cr
 13 \tab Water.Amount    \tab numeric \tab the amount of water added to the cart. \cr
 14 \tab Water.Loss      \tab numeric \tab the amount of water lost since the previous watering. \cr
 15 \tab Area            \tab numeric \tab the total number of plant pixels in three plant images.
}}
\references{
Watts-Williams SJ, Jewell N, Brien C, Berger B, Garnett T, Cavagnaro TR (2019) Using high-throughput phenotyping to explore growth responses to mycorrhizal fungi and zinc in three plant species. \emph{Plant Phenomics}, \bold{2019}, 12.
}
\keyword{datasets}