cat("#### Test corrPlot\n")
test_that("corrPlot", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  library(GGally)
  
  data(exampleData)
  responses <- c("Area","Area.SV","Area.TV", "Image.Biomass", "Max.Height","Centre.Mass",
                 "Density", "Compactness.TV", "Compactness.SV")
 testthat::expect_null(plotCorrmatrix(longi.dat, responses, 
                                        pairs.sets=list(c(1:4),c(5:7))))
})