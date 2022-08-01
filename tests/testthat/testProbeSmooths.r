
cat("#### Test smooth.frames with a baby example\n")
test_that("smooth.frames_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(scales)
dat <- read.table(header = TRUE, text = "
Type TunePar TuneVal Tuning Method       ID  DAP   PSA      sPSA
NCSS      df       4   df-4 direct 045451-C   28 57.446 51.18456
NCSS      df       4   df-4 direct 045451-C   30 89.306 87.67343
NCSS      df       7   df-7 direct 045451-C   28 57.446 57.01589
NCSS      df       7   df-7 direct 045451-C   30 89.306 87.01316
")
  dat[1:7] <- lapply(dat[1:6], factor)
  dat <- as.smooths.frame(dat, individuals = "ID", times = "DAP")
  testthat::expect_true(is.smooths.frame(dat))
  testthat::expect_true(validSmoothsFrame(dat))
  
  tmp <- dat[,-c(2,8)]
  testthat::expect_true(is.smooths.frame(tmp))
  testthat::expect_silent(validsframe <- validSmoothsFrame(tmp))
  testthat::expect_equal(validsframe[2], 
                         paste0("\n  The followng attributes of a smooths.frame are NULL: ",
                                "t, nschemes, individuals, times"))
  testthat::expect_equal(validsframe[3], 
                         "\n  Do not have the following required smoothing-parameters columns in a smooths.frame: TunePar")
  
  testthat::expect_error(as.smooths.frame(tmp), 
                         regexp = paste0("Cannot assign smooths.frame class to supplied data.frame because it does not ",
                                         "contain the following smoothing-parameters columns: TunePar"),                       
                         fixed = TRUE)
  
})


cat("#### Test probeSmooths with Judith Atieno 0278\n")
test_that("chickpea_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  library(ggplot2)
  
  data(dat1)
  
  #'## Values of df for which to obtain plots
  df <- c(4,7)
  
  
  #'## Obtain separate plots Tunings
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             df=df, 
                                             breaks.spacing.x = 2,
                                             plots.by.pf = "Tuning", 
                                             facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse"),
                           regexp = "NaNs produced")
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_equal(ncol(t), 16)
  
  #'## Obtain separate plots Tunings, when includes NCSS for lambda
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             df=df, lambdas = list(NCSS = 0.0001),
                                             plots.by.pf = "Tuning", 
                                             facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse",
                                             include.raw.pf = "alone"),
                           regexp = "NaNs produced")
  testthat::expect_equal(nrow(t), 25440)
  testthat::expect_equal(ncol(t), 16)
  
  
  #'## Obtain separate plots Tunings, when includes NCSS & PS for lambda
  testthat::expect_error(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                           times = "TimeAfterPlanting", 
                                           df=7, lambdas = list(NCSS = 0.0001, PS = 1),
                                           plots.by.pf = "Tuning", 
                                           facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse"),
                         regexp = "The following names for the components of lambdas are not in the specified spline.types: PS")
  
  #'## Obtain separate plots Tunings, when includes NCSS & PS for lambda
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             spline.types = c("NCSS","PS"),
                                             df=7, lambdas = list(NCSS = 0.0001, PS = 1),
                                             plots.by.pf = c("Type","Tuning"), 
                                             facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse",  
                                             include.raw.pf = "facet.x"),
                           regexp = "NaNs produced")
  testthat::expect_equal(nrow(t), 25440)
  testthat::expect_equal(ncol(t), 16)
  
  #'## Obtain separate plots Tunings, when includes NCSS & PS for lambda and set include.raw.pf to "alone"
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             spline.types = c("NCSS","PS"),
                                             df=7, lambdas = list(NCSS = 0.0001, PS = 1),
                                             plots.by.pf = c("Type","Tuning"), 
                                             facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse",  
                                             include.raw.pf = "alone"),
                           regexp = "NaNs produced")
  testthat::expect_equal(nrow(t), 25440)
  testthat::expect_equal(ncol(t), 16)
  
  
  #Test various combinations of smoothing and non-smoothing factors in facet.x.pf
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             get.rates = FALSE, 
                                             keep.columns = c("Smarthouse", "Treatment.1"),
                                             df = df, which.plots = "none"),
                           repexp = "get.rates is FALSE; trait.types changed to response and propn.types.med reduced to its first value")
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_equal(ncol(t), 11)
  testthat::expect_true(all(c("Smarthouse", "Treatment.1") %in% names(t)))
  #test for incorrect trait.types
  testthat::expect_error(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                               times = "TimeAfterPlanting", 
                                               plots.by.pf = "Smarthouse", 
                                               facet.x.pf = c("Treatment.1", "Method", "Tuning")),
                         regexp = paste0("The following traits are not in the smooths.frame: ShootArea1000.AGR, ",
                                         "ShootArea1000.RGR; perhaps, trait.types needs to be set differently"))
  testthat::expect_silent(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                                times = "TimeAfterPlanting", 
                                                trait.types = ("response"), 
                                                plots.by.pf = "Tuning", 
                                                facet.x.pf = c("Smarthouse", "Treatment.1")))
  testthat::expect_silent(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                                times = "TimeAfterPlanting", 
                                                trait.types = ("response"), 
                                                plots.by.pf = "Smarthouse", 
                                                facet.x.pf = c("Method", "Treatment.1", "Tuning")))
  testthat::expect_silent(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                                times = "TimeAfterPlanting", 
                                                trait.types = ("response"), 
                                                plots.by.pf = "Smarthouse", 
                                                facet.x.pf = c("Method", "Treatment.1", "Tuning"), 
                                                collapse.facets.x.pf = FALSE))
  testthat::expect_silent(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                                times = "TimeAfterPlanting", 
                                                trait.types = ("response"), 
                                                plots.by.pf = "Smarthouse", 
                                                facet.x.pf = c("Method", "Treatment.1", "Tuning"),
                                                include.raw.pf = "facet.x"))
  #Don't collapse the facet.x
  testthat::expect_silent(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                                times = "TimeAfterPlanting", 
                                                trait.types = ("response"), 
                                                plots.by.pf = "Smarthouse", 
                                                facet.x.pf = c("Method", "Treatment.1", "Tuning"), 
                                                collapse.facets.x.pf = FALSE))
  #Inclusion of raw with both smoothing-parameter and other factors in facet.x
  testthat::expect_silent(
    plotSmoothsComparison(data = t, response="ShootArea1000", 
                          times = "TimeAfterPlanting", 
                          trait.types = ("response"), 
                          plots.by.pf = "Smarthouse", 
                          facet.x.pf = c("Method", "Treatment.1", "Tuning"), 
                          collapse.facets.x.pf = FALSE, include.raw.pf = "facet.x"))
  #Inclusion of raw with both smoothing-parameter and other factors in facet.y
  testthat::expect_silent(
    plotSmoothsComparison(data = t, response="ShootArea1000", 
                          times = "TimeAfterPlanting", 
                          trait.types = ("response"), 
                          plots.by.pf = "Smarthouse", 
                          facet.y.pf = c("Method", "Treatment.1", "Tuning"), 
                          collapse.facets.y.pf = FALSE, include.raw.pf = "facet.y"))
  #Inclusion of raw with other factors in facet.x and  facet.y
  testthat::expect_silent(
    plotSmoothsComparison(data = t, response="ShootArea1000", 
                          times = "TimeAfterPlanting", 
                          trait.types = ("response"), 
                          plots.by.pf = "Smarthouse", 
                          facet.x.pf = c("Method", "Tuning"), facet.y.pf = "Treatment.1",  
                          collapse.facets.x.pf = FALSE, include.raw.pf = "facet.x"))
  #removes Median whiskers
  testthat::expect_silent(
    plotSmoothsComparison(data = t, response="ShootArea1000", 
                          times = "TimeAfterPlanting", 
                          trait.types = ("response"), 
                          plots.by.pf = "Smarthouse", 
                          facet.y.pf = c("Method", "Treatment.1", "Tuning"), 
                          collapse.facets.y.pf = FALSE, include.raw.pf = "facet.y"))
   #removes Median whiskers
  testthat::expect_silent(plotSmoothsComparison(data = t, response="ShootArea1000", 
                                                times = "TimeAfterPlanting", 
                                                trait.types = ("response"), 
                                                plots.by.pf = "Smarthouse", 
                                                facet.x.pf = c("Method", "Treatment.1", "Tuning"),
                                                addMediansWhiskers = FALSE))
  #Must specify either plots.by.pf or plots.compare 
  testthat::expect_warning(testthat::expect_error(
    t <- probeSmooths(dat1, response="ShootArea1000", 
                      times = "TimeAfterPlanting", df=df),
    regexp = paste0("There are no smoothing-parameter factors assigned to the plots and facet ", 
                    "arguments - enough of them need to be assigned so that they uniquely index ", 
                    "the combinations of the smoothing-parameter values in the smooths.frame")),
    regexp = "NaNs produced")
  
  #Plot medians.deviations only with only facet.y.med set
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             trait.types = "response", 
                                             df=df, 
                                             which.plots = "medians.dev", 
                                             facet.x.pf = "Tuning", 
                                             facet.y.pf = c("Smarthouse", "Treatment.1"), 
                                             propn.types.med = 0.05,
                                             plots.group.med = "Tuning", 
                                             facet.y.med = c("Smarthouse", "Treatment.1")))
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_equal(ncol(t), 11)
  
  #Plot with plots.by.pf set to Tuning, plots.compare set to NULL and facet.x.pf set to Treatment.1
  testthat::expect_silent(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                            times = "TimeAfterPlanting", 
                                            trait.types=c("AGR"),
                                            df=df, 
                                            plots.by.pf = "Tuning", 
                                            facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse", 
                                            include.raw.pf = "alone"))
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_equal(ncol(t), 14)
  
  #specify Treatment.1 and Tuning for facet.x.pf
  testthat::expect_silent(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                            times = "TimeAfterPlanting", 
                                            trait.types=c("response"), get.rates=FALSE, 
                                            df=df, 
                                            facet.x.pf = c("Treatment.1", "Tuning"), 
                                            facet.y.pf = "Smarthouse"))
  
  #include raw plots with compare and no facet.x.pf
  testthat::expect_silent(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                            times = "TimeAfterPlanting", 
                                            trait.types=c("response"), get.rates=FALSE, 
                                            df=df, 
                                            facet.x.pf = "Tuning", facet.y.pf = "Treatment.1", 
                                            include.raw.pf = "facet.x"))
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_true(all(c("Type", "TunePar", "TuneVal", "Tuning", "Method", 
                              "Snapshot.ID.Tag", "TimeAfterPlanting",  
                              "Treatment.1", "ShootArea1000", "sShootArea1000") 
                            %in% names(t)))
  
  
  
  #'## Plots of response and deviations boxplots 
  testthat::expect_silent(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                            times = "TimeAfterPlanting", 
                                            x.title = "DAP", 
                                            trait.types=c("response"), get.rates=FALSE, 
                                            which.plots = c("profiles", "absolute"),
                                            df=df, 
                                            facet.x.pf = "Tuning", 
                                            facet.y.pf = "Treatment.1"))
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_true(all(c("Type", "TunePar", "TuneVal", "Tuning", "Method", 
                              "Snapshot.ID.Tag", "TimeAfterPlanting",  
                              "Treatment.1", "ShootArea1000", "sShootArea1000") 
                            %in% names(t)))
  
  #plots.by.pf gives separate plots
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             trait.types=c("response", "AGR"), 
                                             x.title = "DAP", 
                                             df=df, 
                                             plots.by.pf = "Tuning", 
                                             facet.x.pf = "Treatment.1", facet.y.pf = "Smarthouse", 
                                             which.plots = c("profiles", "absolute", "relative")))
  testthat::expect_equal(nrow(t), 16960)
  testthat::expect_lt(max(abs(t[t$TuneVal == "4", "ShootArea1000"] - t[t$TuneVal == "4", "sShootArea1000"]), 
                          na.rm = TRUE) - 5.563407, 1e-07)
  testthat::expect_lt(max(abs(t[t$TuneVal == "4", "ShootArea1000"] - t[t$TuneVal == "4", "sShootArea1000"])/
                            t[t$TuneVal == "4", "sShootArea1000"]) - 370.7951, 1e-04)
  
  testthat::expect_warning(t <- probeSmooths(data = dat1, response="ShootArea1000", 
                                             times = "TimeAfterPlanting", 
                                             x.title = "DAP", 
                                             trait.types=c("response", "AGR"), 
                                             df=df, 
                                             plots.by.pf = "Tuning", 
                                             facet.y.pf = c("Smarthouse", "Treatment.1"),
                                             which.plots = c("profiles", "absolute", "relative")))
  testthat::expect_equal(nrow(t), 16960)
  
})

cat("#### Test probeSmooths with small example\n")
test_that("exampleData_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(growthPheno)
  
  data(exampleData)
  
  vline <- list(ggplot2::geom_vline(xintercept=29, linetype="longdash", size=1))
  
  plotDeviationsBoxes(longi.dat, observed = "PSA", smoothed = "sPSA",
                      x.factor="DAP", df =5)
  
  testthat::expect_silent(tmp <- probeSmooths(data = longi.dat, 
                                              response = "PSA", times = "DAP", 
                                              df = c(4,7), 
                                              facet.x.pf = "Tuning", facet.y.pf = "Treatment.1",
                                              ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(tmp), 560)
  testthat::expect_equal(ncol(tmp), 15)
  
  testthat::expect_warning(testthat::expect_error(
    tmp <- probeSmooths(data = longi.dat, 
                        response = "PSA", times = "DAP", 
                        df = c(4,7), 
                        facet.x.pf = "Method", 
                        get.rates = FALSE,
                        ggplotFuncsProfile = vline),
    regexp = paste0("\n The number of different combinations of the smoothing-parameter values in the smooths.frame and ", 
                    "of the levels combination of the following factors nominated in the facet/plots arguments ", 
                    "are not equal: Method")),
    regexp = paste0("get.rates is FALSE; trait.types changed to response and propn.types.med reduced to its first value"))
  testthat::expect_warning(tmp <- probeSmooths(data = longi.dat, 
                                               response = "PSA", times = "DAP", 
                                               df = c(4,7), 
                                               facet.x.pf = "Tuning", 
                                               get.rates = FALSE,
                                               ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(tmp), 560)
  testthat::expect_equal(ncol(tmp), 9)
  
  testthat::expect_warning(tmp <- probeSmooths(data = longi.dat, 
                                               response = "PSA", times = "DAP", 
                                               df = c(4:7), 
                                               facet.x.pf = "Tuning", 
                                               alpha.pf = 0.6, get.rates = FALSE,
                                               ggplotFuncsProfile = vline),
                           regexp = paste0("get.rates is FALSE; trait.types changed to response and ",
                                           "propn.types.med reduced to its first value"))
  testthat::expect_equal(nrow(tmp), 1120)
  testthat::expect_equal(ncol(tmp), 9)
  
  #AGR trait plot only
  testthat::expect_silent(tmp <- probeSmooths(data = longi.dat, 
                                              response = "PSA", times = "DAP", 
                                              df = c(4:7), 
                                              facet.x.pf = "Tuning", 
                                              alpha.pf = 0.5, trait.types = "AGR",
                                              ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(tmp), 1120)
  testthat::expect_equal(ncol(tmp), 12)
  
  #Compare medians and profiles for multiple traits  
  testthat::expect_warning(tmp <- probeSmooths(data = longi.dat, 
                                               response = "PSA", times = "DAP", 
                                               df = c(4:7), 
                                               facet.x.pf = "Tuning", 
                                               which.plots = c("profiles", "medians.dev"), 
                                               plots.group.med = "Tuning", 
                                               propn.types.med = c(0.02,0.1, 0.2),
                                               ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(tmp), 1120)
  testthat::expect_equal(ncol(tmp), 14)
  
  #Compare medians for multiple traits  
  testthat::expect_warning(tmp <- probeSmooths(data = longi.dat, 
                                               response = "PSA", times = "DAP", 
                                               df = c(4:7), 
                                               which.plots = "medians.dev", 
                                               plots.group.med = "Tuning",
                                               propn.types.med = c(0.02,0.1, 0.2),
                                               ggplotFuncsProfile = vline))
  testthat::expect_equal(nrow(tmp), 1120)
  testthat::expect_equal(ncol(tmp), 14)
  
  testthat::expect_silent(tmp <- probeSmooths(data = longi.dat, 
                                              response = "PSA", times = "DAP", 
                                              df = c(4,7), 
                                              plots.by.pf = "Tuning", facet.y.pf = "Treatment.1",
                                              ggplotFuncsProfile = vline))
  
  #Test for a single line per plot - caused by plots.by.med
  testthat::expect_silent(med <- plotSmoothsMedianDevns(data = tmp, 
                                                        response = "PSA", 
                                                        response.smoothed = "sPSA", 
                                                        times = "DAP", 
                                                        trait.types = "response", 
                                                        x.title = "DAP", 
                                                        plots.by.med = "Tuning",  
                                                        propn.types.med = 0.05,
                                                        ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med$plots), 1)
  testthat::expect_equal(nrow(med$med.devn.dat), 28)
  testthat::expect_equal(ncol(med$med.devn.dat), 3)

  #Test for a single line per plot - caused by plots.by.med; specification of colour
  testthat::expect_silent(med <- plotSmoothsMedianDevns(data = tmp, 
                                                        response = "PSA", 
                                                        response.smoothed = "sPSA", 
                                                        times = "DAP", 
                                                        trait.types = "response", 
                                                        x.title = "DAP", 
                                                        plots.by.med = "Tuning", 
                                                        colour.values.med = "blue", shape.values.med = 17,
                                                        propn.types.med = 0.05,
                                                        ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med$plots), 1)
  testthat::expect_equal(nrow(med$med.devn.dat), 28)
  testthat::expect_equal(ncol(med$med.devn.dat), 3)
  
  #Compare medians and absolute deviations for multiple traits  
  testthat::expect_warning(traits <- probeSmooths(data = longi.dat, 
                                                  response = "PSA", times = "DAP", 
                                                  df = c(4:7), 
                                                  which.plots = c("medians.dev", "absolute"),
                                                  facet.x.pf = "Tuning", #for the absolute.deviations
                                                  plots.group.med = "Tuning",  
                                                  propn.types.med = NULL))
  testthat::expect_equal(nrow(traits), 1120)
  testthat::expect_equal(ncol(traits), 14)
  
  #Form and save the data.frame containing the smooths produced by probeSmooths for testing associated functions
  testthat::expect_silent(traits <- probeSmooths(data = longi.dat, 
                                                 response = "PSA", times = "DAP", 
                                                 df = c(4:7), 
                                                 which.plots = "none", 
                                                 keep.columns = c("Treatment.1", "Genotype.ID"))) #so in smooths.frame
  
  #Fit simple and four-parameter logistics and obtain the fitted values
  library(nlme)
  extra.dat <- longi.dat[, -grep("sPSA", names(longi.dat), fixed = TRUE)]  
  logist.grp <- nlme::groupedData(PSA ~ xDAP | Snapshot.ID.Tag, data = longi.dat)
  #Fit the simple logistic model  
  logist.lis <- nlme::nlsList(SSlogis, logist.grp, na.action = na.pass)
  logist.dat <- within(extra.dat, sPSA <- fitted(logist.lis))
  logist.dat <- byIndv4Times_GRsDiff(data = logist.dat, 
                                     response = "sPSA", 
                                     individuals="Snapshot.ID.Tag",
                                     times="DAP", 
                                     which.rates = c("AGR", "RGR"))
  logist.dat <- cbind(Tuning = factor("Logistic"), logist.dat)
  #Fit the four-parameter logistic model - generates warnings
  logis4.lis <- suppressWarnings(nlme::nlsList(SSfpl, logist.grp, na.action = na.pass))
  logis4.dat <- within(extra.dat, sPSA <- fitted(logis4.lis))
  logis4.dat <- byIndv4Times_GRsDiff(data = logis4.dat, 
                                     response = "sPSA", 
                                     individuals="Snapshot.ID.Tag",
                                     times="DAP", 
                                     which.rates = c("AGR", "RGR"))
  logis4.dat <- cbind(Tuning = factor("Logis-4par"), logis4.dat)
  #Combine the logistic fits
  extra.dat <- rbind(logist.dat,logis4.dat)
  extra.dat <- cbind(Type = factor("NonLinear"), extra.dat)
  
 
  #Check computation of median deviations for Control Treatment, Genotype = 120855 and df = 5
  tmp <- subset(traits, Treatment.1 == "Control" & Genotype.ID == "120855" & TuneVal == "4")
  tmp$PSA.devn <- tmp$PSA - tmp$sPSA
  med.vals <- tapply(tmp$PSA.devn, tmp$DAP, median, na.rm = TRUE)
  testthat::expect_silent(med <- plotSmoothsMedianDevns(data = traits, 
                                                        response = "PSA", times = "DAP", 
                                                        x.title = "DAP", 
                                                        trait.types = "response", 
                                                        plots.group.med = "Tuning",  
                                                        facet.y.med = c("Treatment.1", "Genotype.ID"),
                                                        propn.types.med = 0.05,
                                                        ggplotFuncsMedDevn = vline, printPlot = FALSE))
  med.dat <- med$med.devn.dat
  med.dat <- subset(med.dat, Treatment.1 == "Control" & Genotype.ID == "120855" & SmoothParams == "df-4")
  testthat::expect_true(all(abs(med.dat$PSA.devn - med.vals) < 1e-5))
  
  testthat::expect_silent(med <- plotSmoothsMedianDevns(data = traits, 
                                                        response = "PSA", 
                                                        response.smoothed = "sPSA", 
                                                        times = "DAP", 
                                                        trait.types = "response", 
                                                        x.title = "DAP", 
                                                        plots.group.med = "Tuning",  
                                                        propn.types.med = 0.05,
                                                        ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(nrow(med$med.devn.dat), 56)
  testthat::expect_equal(ncol(med$med.devn.dat), 3)
  testthat::expect_equal(length(med$plots), 1)
  
  
  #Test external.smooths with probeSmooths
  testthat::expect_warning(smth.extra <- probeSmooths(data = longi.dat, 
                                                      response = "PSA", times = "DAP", 
                                                      df = c(4:7), 
                                                      facet.x.pf = "Tuning", 
                                                      which.plots = c("medians.dev", "absolute"),
                                                      plots.group.med = "Tuning",
                                                      propn.types.med = NULL,
                                                      external.smooths = extra.dat))
  testthat::expect_equal(nrow(smth.extra), 1680)
  testthat::expect_equal(ncol(smth.extra), 14)
  
  #Use external.smooths incorporated using probeSmooths
  testthat::expect_warning(med <- plotSmoothsMedianDevns(data = smth.extra, 
                                                         response = "PSA", times = "DAP", 
                                                         x.title = "DAP", 
                                                         plots.by.med = "Type", 
                                                         plots.group.med = "Tuning",  
                                                         propn.types.med = c(0.02,0.1, 0.2),
                                                         ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(nrow(med$med.devn.dat), 84)
  testthat::expect_equal(ncol(med$med.devn.dat), 6)
  testthat::expect_equal(levels(med$med.devn.dat$SmoothParams), 
                         c("df-4","df-5","df-6","df-7","Logistic","Logis-4par"))
  testthat::expect_equal(length(med$plots), 3)
  
  #Form a smooths.frame with external.smooths for further testing
  testthat::expect_silent(smth.extra <- probeSmooths(data = longi.dat, 
                                                     response = "PSA", times = "DAP", 
                                                     df = c(4:7), 
                                                     facet.x.pf = "Tuning", 
                                                     include.raw.pf = "facet.x",
                                                     propn.types.med = NULL,
                                                     external.smooths = extra.dat))
  testthat::expect_equal(nrow(smth.extra), 1680)
  testthat::expect_equal(ncol(smth.extra), 14)
  #Use plots.by.pf different to that in  probeSmooths
  plts <- plotSmoothsComparison(data = smth.extra, 
                                response = "PSA", times = "DAP", 
                                x.title = "DAP", include.raw.pf = "facet.x", 
                                plots.by.pf = "Type", facet.x.pf = "Tuning",
                                alpha.pf = 0.2, 
                                ggplotFuncsProfile = vline)  
  testthat::expect_equal(length(plts), 3)
  testthat::expect_true(all(c("PSA","PSA.AGR","PSA.RGR") %in% names(plts)))
  testthat::expect_true(all(c("profiles", "deviations") %in% names(plts$PSA)))
  testthat::expect_true(all(c("profiles", "deviations") %in% names(plts$PSA.AGR)))
  testthat::expect_equal(length(plts$PSA$deviations), 0)
  testthat::expect_equal(length(plts$PSA$profiles), 2)
  testthat::expect_true(all(c("NCSS", "NonLinear") %in% names(plts$PSA$profiles)))
  
  #Test addMediansWhiskers
  plotSmoothsComparison(data = smth.extra, response = "PSA", times = "DAP", 
                        x.title = "DAP", include.raw.pf = "facet.x", 
                        facet.x.pf = c("Type", "Tuning"),
                        colour.column.pf = "Method", 
                        colour.values.pf = c("olivedrab", "orange"),
                        alpha.pf = 0.25, addMediansWhiskers = TRUE,
                        ggplotFuncsProfile = vline)  
  
  plts <- plotSmoothsComparison(data = smth.extra, response = "PSA", times = "DAP", 
                                x.title = "DAP", include.raw.pf = "facet.x", 
                                plots.by.pf = "Type", facet.x.pf = "Tuning",
                                colour.pf = "olivedrab",
                                alpha.pf = 0.25, addMediansWhiskers = TRUE,
                                ggplotFuncsProfile = vline)  
  
  
  #Add external smooths and include deviations plots
  testthat::expect_warning(smth.extra <- probeSmooths(data = longi.dat, 
                                                      response = "PSA", times ="DAP", 
                                                      df = c(4:7), 
                                                      facet.x.pf = c("Type", "Tuning"),  
                                                      which.plots = c("medians.dev", "absolute"),
                                                      plots.group.med = c("Type", "Tuning"),  
                                                      propn.types.med = NULL,
                                                      external.smooths = extra.dat))
  testthat::expect_equal(nrow(smth.extra), 1680)
  testthat::expect_equal(ncol(smth.extra), 14)
  
  
  testthat::expect_warning(med <- plotSmoothsMedianDevns(data = smth.extra, 
                                                         response = "PSA", times = "DAP", 
                                                         x.title = "DAP", 
                                                         plots.group.med = c("Type", "Tuning"),  
                                                         propn.types.med = c(0.02,0.1, 0.2),
                                                         ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(nrow(med$med.devn.dat), 84)
  testthat::expect_equal(ncol(med$med.devn.dat), 5)
  testthat::expect_equal(levels(med$med.devn.dat$SmoothParams), 
                         c("NCSS-df-4","NCSS-df-5","NCSS-df-6","NCSS-df-7",
                           "NonLinear-Logistic","NonLinear-Logis-4par"))
  testthat::expect_equal(length(med$plots), 3)
  
  testthat::expect_warning(med <- plotSmoothsMedianDevns(data = smth.extra, 
                                                         response = "PSA", times = "DAP", 
                                                         x.title = "DAP", 
                                                         facet.x.med = "Type",
                                                         plots.group.med = "Tuning",  
                                                         propn.types.med = c(0.02,0.1, 0.2),
                                                         ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(levels(med$med.devn.dat$SmoothParams), 
                         c("df-4","df-5","df-6","df-7","Logistic","Logis-4par"))
  
  testthat::expect_warning(med <- plotSmoothsMedianDevns(data = smth.extra, 
                                                         response = "PSA", times = "DAP", 
                                                         x.title = "DAP", 
                                                         plots.group.med = "Tuning",  
                                                         facet.y.med = "Type",
                                                         propn.types.med = c(0.02,0.1, 0.2),
                                                         ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(levels(med$med.devn.dat$SmoothParams), 
                         c("df-4","df-5","df-6","df-7","Logistic","Logis-4par"))
  
  testthat::expect_warning(med <- plotSmoothsMedianDevns(data = smth.extra, 
                                                         response = "PSA", times = "DAP", 
                                                         x.title = "DAP", 
                                                         plots.group.med = "Tuning",  
                                                         facet.x.med = "Type",
                                                         propn.types.med = c(0.02,0.1, 0.2),
                                                         ggplotFuncsMedDevn = vline))
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(levels(med$med.devn.dat$SmoothParams), 
                         c("df-4","df-5","df-6","df-7","Logistic","Logis-4par"))
})


cat("#### Test probeSmooths with tomato example\n")
test_that("tomato_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(growthPheno)
  
  data(tomato.dat)

  df.vec        <- c(4:6,12)
  labelMyc <- as_labeller(function(lev) paste(lev, "AMF"))
  labelZn <- as_labeller(function(lev) paste("Zn:", lev, "ppm"))
  DAP.endpts   <- c(18,22,27,33,39,43,51)
  nDAP.endpts <- length(DAP.endpts)
  DAP.starts <- DAP.endpts[-nDAP.endpts]
  DAP.stops   <- DAP.endpts[-1]
  
  #'## Gives error that the Length of propn.types.med is not the same as the number of traits
  testthat::expect_warning(
    tom <- probeSmooths(data = tomato.dat, 
                        response = "PSA", response.smoothed = "sPSA", 
                        times = "DAP", 
                        smoothing.methods = c("dir", "log"), 
                        df = c(4,7), get.rates = FALSE,
                        which.plots = c("profiles", "medians.dev"),
                        facet.y.pf = c("Zn","AMF"),
                        plots.by.pf = "Tuning", facet.x.pf = "Method",
                        plots.by.med = "Tuning", plots.group.med = "Method", 
                        facet.y.med = c("Zn","AMF"),
                        labeller = labeller(Zn = labelZn, 
                                            AMF = labelMyc)),
    regexp = "get.rates is FALSE; trait.types changed to response and propn.types.med reduced to its first value")
  testthat::expect_equal(nrow(tom), 4480)
  testthat::expect_equal(ncol(tom), 11)
  
  testthat::expect_silent(med <- plotSmoothsMedianDevns(data = tom, 
                                                        response = "PSA", 
                                                        response.smoothed = "sPSA", 
                                                        times = "DAP", 
                                                        trait.types = "response", 
                                                        x.title = "DAP", 
                                                        y.titles = "PSA deviation (kpixels)",
                                                        facet.x.med = "Zn", 
                                                        facet.y.med = "AMF",
                                                        plots.group.med = c("Method","Tuning"), 
                                                        propn.types.med = 0.1,
                                                        labeller = labeller(Zn = labelZn, 
                                                                            AMF = labelMyc)))
  
  testthat::expect_equal(length(med$plots), 1)
  testthat::expect_equal(nrow(med$med.devn.dat), 1120)
  testthat::expect_equal(ncol(med$med.devn.dat), 5)
  
  #Multiple df, single methods
  testthat::expect_warning(tom <- probeSmooths(data = tomato.dat, response = "PSA", 
                                               response.smoothed = "sPSA", 
                                               which.plots = "medians.dev", 
                                               times = "DAP",  df=5:6,
                                               smoothing.methods = c("logarithmic"),
                                               propn.types.med = c(0.02, 0.2, 0.5),
                                               plots.group.med = "Tuning", plots.by.med = "Method"))
  testthat::expect_equal(nrow(tom), 2240)
  testthat::expect_equal(ncol(tom), 14)
  
  #'Single `df`, multiple methods and trait.types
  testthat::expect_warning(
    tomdiff <- probeSmooths(data = tomato.dat, response = "PSA", 
                            response.smoothed = "sPSA", 
                            times = "DAP", df=5,
                            smoothing.methods = c("direct","logarithmic"),
                            which.plots = "none"))
  testthat::expect_equal(nrow(tomdiff), 2240)
  testthat::expect_equal(ncol(tomdiff), 14)
  testthat::expect_warning(med <- plotSmoothsMedianDevns(data = tomdiff, response = "PSA", 
                                                         response.smoothed = "sPSA", 
                                                         times = "DAP", x.title = "DAP", 
                                                         plots.by.med = "Tuning", 
                                                         plots.group.med = "Method", 
                                                         propn.types.med = c(0.02, 0.2, 0.5)))
  
  testthat::expect_equal(length(med), 2)
  testthat::expect_true(all(names(med) == c("plots", "med.devn.dat")))
  testthat::expect_equal(nrow(med$med.devn.dat), 70)
  testthat::expect_equal(ncol(med$med.devn.dat), 6)
  testthat::expect_equal(length(med$plots), 3)
  testthat::expect_true(all(names(med$plots) == c("PSA","PSA.AGR","PSA.RGR")))
  testthat::expect_warning(print(med$plots$PSA.AGR$`df-5`))
  
  #'Single `df`, single method - plots.by.med = Tuning
  testthat::expect_warning(tom <- probeSmooths(data = tomato.dat, response = "PSA", 
                                               response.smoothed = "sPSA", 
                                               times = "DAP", df=5,
                                               smoothing.methods = c("direct"),
                                               which.plots = "medians.dev", 
                                               propn.types.med = c(0.02, 0.2, 0.5),
                                               plots.by.med = "Tuning",
                                               colour.values.med = "blue"))
  testthat::expect_equal(nrow(tom), 1120)
  testthat::expect_equal(ncol(tom), 14)
  
  #'Single `df`, single method - plots.group.med = Tuning
  testthat::expect_warning(tom <- probeSmooths(data = tomato.dat, response = "PSA", 
                                               response.smoothed = "sPSA", 
                                               times = "DAP", df=5,
                                               smoothing.methods = c("direct"),
                                               which.plots = "medians.dev", 
                                               propn.types.med = c(0.02, 0.2, 0.5),
                                               plots.group.med = "Tuning"))
  testthat::expect_equal(nrow(tom), 1120)
  testthat::expect_equal(ncol(tom), 14)
  
  #test custom schemes
  spar.schemes <- data.frame(Type = c("N", "NCS", "P"), 
                             TunePar = c("df", "df", "lam"),
                             TuneVal = c(4, 6, 1),
                             Method = c("dir", "log", "log"))
  cols <- scales::brewer_pal("div", "Paired")(6)[c(2,4,6,8,1,3,5,7)]
  testthat::expect_warning(tom <- probeSmooths(data = tomato.dat, response = "PSA", 
                                               response.smoothed = "sPSA", 
                                               times = "DAP", 
                                               smoothing.schemes = spar.schemes,
                                               which.plots = "medians.dev", 
                                               colour.values.med = cols[1:3], 
                                               propn.types.med = c(0.02, 0.2, 0.5),
                                               plots.group.med = c("Type", "Tuning", "Method")))
  testthat::expect_equal(nrow(tom), 3360)
  testthat::expect_equal(ncol(tom), 14)
  
  #Multiple df, two spline.types - uses npspline.segments
  testthat::expect_warning(tom <- probeSmooths(data = tomato.dat, response = "PSA", 
                                               response.smoothed = "sPSA", 
                                               times = "DAP",  
                                               which.plots = c("medians.dev", "profiles"), 
                                               smoothing.methods = "log", spline.types = c("NCSS", "PS"), 
                                               df=5:6, lambdas = c(0.1,1), npspline.segments = 30,
                                               facet.x.pf = c("Type", "Tuning"),
                                               propn.types.med = c(0.02, 0.2, 0.5),
                                               plots.group.med = "Tuning", facet.x.med = "Type"))
  testthat::expect_equal(nrow(tom), 4480)
  testthat::expect_equal(ncol(tom), 14)
  
  
  
  #probeSmooths without deviations plots for next set of tests
  testthat::expect_warning(tom <- probeSmooths(data = tomato.dat, response = "PSA", 
                                               response.smoothed = "sPSA", 
                                               times = "DAP", 
                                               smoothing.schemes = spar.schemes,
                                               which.plots = "none"),
                           regexp = "NaNs produced")
  
  #test various combinations of plots.by.pf, plots.compare and include.raw.pf alone
  testthat::expect_silent(plts <- 
                            plotSmoothsComparison(data = tom, response = "PSA", 
                                                  response.smoothed = "sPSA", 
                                                  times = "DAP", 
                                                  which.plots = "profiles", 
                                                  include.raw.pf = "alone",
                                                  plots.by.pf = c("Type", "Tuning", "Method"),
                                                  colour.column.pf = "Method", 
                                                  colour.values.pf = c("orange", "olivedrab"), 
                                                  addMediansWhiskers = TRUE))
  testthat::expect_equal(length(plts$PSA$profiles), 4)
  testthat::expect_equal(length(plts$PSA$deviations), 0)
  

  #test various combinations of plots.by.pf, plots.compare and include.raw.pf to facet.x without facet.x.pf set
  testthat::expect_error(
    plotSmoothsComparison(data = tom, response = "PSA", 
                          response.smoothed = "sPSA", 
                          times = "DAP", 
                          which.plots = "profiles", 
                          include.raw.pf = "facet.x",
                          plots.by.pf = c("Type", "Tuning", "Method"),
                          colour.column.pf = "Method", 
                          colour.values.pf = c("orange", "olivedrab"), 
                          addMediansWhiskers = TRUE), 
    regexp = "The argument incl.raw.pf is set to facet.x, but facet.x.pf has not been set to include a variable")
  testthat::expect_equal(length(plts$PSA$profiles), 4)
  testthat::expect_equal(length(plts$PSA$deviations), 0)
  
  testthat::expect_silent(plotSmoothsComparison(data = tom, response = "PSA", 
                                                response.smoothed = "sPSA", 
                                                times = "DAP", 
                                                which.plots = "profiles", 
                                                include.raw.pf = "no",
                                                plots.by.pf = c("Type", "Tuning", "Method")))
  
  testthat::expect_silent(plotSmoothsComparison(data = tom, response = "PSA", 
                                                response.smoothed = "sPSA", 
                                                times = "DAP", 
                                                which.plots = "profiles", 
                                                include.raw.pf = "facet.x",
                                                facet.x.pf = c("Type","Tuning","Method"),
                                                alpha.pf = 0.4, colour.column.pf = "Method", 
                                                colour.values.pf = c("orange", "olivedrab"), 
                                                addMediansWhiskers = TRUE))
  
  testthat::expect_silent(plotSmoothsComparison(data = tom, response = "PSA", 
                                                response.smoothed = "sPSA", 
                                                times = "DAP", 
                                                which.plots = "profiles", 
                                                include.raw.pf = "facet.x",
                                                colour.pf = "orange", 
                                                facet.x.pf = "Tuning"))
  
  testthat::expect_silent(plotSmoothsComparison(data = tom, response = "PSA", 
                                                response.smoothed = "sPSA", 
                                                times = "DAP", 
                                                which.plots = "profiles", 
                                                include.raw.pf = "no",
                                                facet.x.pf = c("Type","Tuning","Method")))
  
  #'Single `df`, multiple methods and trait.types and GRs using deriv
  suppressWarnings(
    tomdv <- probeSmooths(data = tomato.dat, response = "PSA", 
                          response.smoothed = "sPSA", 
                          times = "DAP", 
                          df=5,smoothing.methods = c("direct","logarithmic"), 
                          rates.method = "deriv",
                          which.plots = "profiles",
                          facet.x.pf = "Method", include.raw.pf = "facet.x"))
  testthat::expect_equal(nrow(tomdv), 2240)
  testthat::expect_equal(ncol(tomdv), 14)
  testthat::expect_true(all(abs(tomdiff$sPSA-tomdv$sPSA) < 1e-08))
  
  tomspl <- byIndv4Times_SplinesGRs(data = tomato.dat, 
                                    response = "PSA", response.smoothed = "sPSA", 
                                    times="DAP", 
                                    df = 5, rates.method = "deriv", 
                                    which.rates = c("AGR", "RGR"), 
                                    suffices.rates = c("AGRdv", "RGRdv"))
  testthat::expect_true(all(abs(tomspl$sPSA-tomdv$sPSA[1:1120]) < 1e-08))
  testthat::expect_true(all(abs(tomspl$sPSA.AGRdv-tomdv$sPSA.AGR[1:1120]) < 1e-08))
  testthat::expect_true(all(abs(tomdv$sPSA.AGR[1:4] - 
                                  c(4.483901, 4.533483, 4.706856, 5.038068)) < 1e-05))
  testthat::expect_true(all(abs(tomspl$sPSA.RGRdv-tomdv$sPSA.RGR[1:1120]) < 1e-08))
  testthat::expect_true(all(abs(tomdv$sPSA.RGR[1:4] - 
                                  c(17.6807353, 0.9536113, 0.5027707, 0.3542860)) < 1e-05))
  
  #Look at water.use
  lambdas <- round(10^c(-0.5, 0, 0.5, 1), digits = 3)
  df = c(4:6)
  x.axis <- list(theme(axis.text.x = element_text(angle = 90),
                       panel.grid.minor.x = element_blank()))
  vline.DAP.endpts <- list(geom_vline(xintercept=DAP.starts, linetype="longdash", 
                                      alpha = 0.5, size=0.75))
  theme.profile <- list(vline.DAP.endpts,x.axis)
  tom.H2O <- probeSmooths(data = tomato.dat, response = "WU", 
                          response.smoothed = "sWU", 
                          times = "DAP",  get.rates = FALSE, 
                          smoothing.methods = "dir", spline.types = c("NCSS", "PS"), 
                          df=df, lambdas = list(PS = lambdas), 
                          which.plots = c("medians.dev", "profiles"), 
                          include.raw.pf = "facet.x",
                          plots.by.pf = "Type", 
                          facet.x.pf = c("Tuning"),
                          propn.types.med = NULL,
                          plots.group.med = "Tuning", facet.x.med = "Type",
                          ggplotFuncsProfile = theme.profile)
  testthat::expect_equal(nrow(tom.H2O), 7840)
  testthat::expect_equal(ncol(tom.H2O), 9)
  
  
})



cat("#### Test probeSmooths Rice experiment\n")
test_that("RicePrepped_growthPheno", {
  skip_if_not_installed("growthPheno")
  skip_on_cran()
  library(dae)
  library(growthPheno)
  
  data(RicePrepped.dat)
  
  testthat::expect_warning(
    smth <- probeSmooths(data = RicePrepped.dat, response = "PSA", response.smoothed = "sPSA", 
                         times = "DAST", 
                         spline.types = c("NCSS", "PS"), smoothing.methods = "log",
                         df = c(4,7), lambdas = c(0.1,1,10), 
                         keep.columns = c("Smarthouse","Salinity"), #so these are included in the smooths.frame
                         which.plots = "none"),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(nrow(smth), 73920)
  testthat::expect_equal(ncol(smth), 16)
  
  #Check that PS smoothing is correct 
  oneplant <- subset(smth, Snapshot.ID.Tag == "045451-C" & Type == "PS" & TuneVal == "0.1")
  fity <- JOPS::psNormal(x = oneplant$DAST, y = log(oneplant$PSA), nseg = 10, lambda = 0.1, 
                         xgrid = oneplant$DAST)
  testthat::expect_true(all(abs(oneplant$sPSA-exp(fity$muhat[,1])) < 1e-05))
  
  #Tuning only for PSA
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.group.med = "Tuning", 
                                facet.y.med = c("Smarthouse","Salinity"),
                                propn.types.med = 0.025))
  testthat::expect_equal(length(t$plots), 1)
  testthat::expect_equal(nrow(t$med.devn.dat), 280)
  testthat::expect_equal(ncol(t$med.devn.dat), 5)
  testthat::expect_true(all(c("Smarthouse", "Salinity", "SmoothParams", "PSA.devn", "DAST") 
                            %in% names(t$med.devn.dat)))
  #Use plots.by.med for Type
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.by.med = "Type", plots.group.med = "Tuning", 
                                facet.y.med = c("Smarthouse","Salinity"),
                                propn.types.med = 0.025))
  testthat::expect_equal(length(t$plots), 1)
  testthat::expect_equal(nrow(t$med.devn.dat), 280)
  testthat::expect_equal(ncol(t$med.devn.dat), 6)
  testthat::expect_true(all(c("fac.by","Smarthouse", "Salinity", "SmoothParams", "PSA.devn", "DAST") 
                            %in% names(t$med.devn.dat)))
  #Use facet.x.med for Type, additional to plots.group 
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.group.med = "Tuning", 
                                facet.x.med = "Type", 
                                facet.y.med = c("Smarthouse","Salinity"),
                                propn.types.med = 0.025))
  testthat::expect_equal(length(t$plots), 1)
  testthat::expect_equal(nrow(t$med.devn.dat), 280)
  testthat::expect_equal(ncol(t$med.devn.dat), 6)
  testthat::expect_true(all(c("Type","Smarthouse", "Salinity", "SmoothParams", "PSA.devn", "DAST") 
                            %in% names(t$med.devn.dat)))
  
  #facet.x a mixture of a smoothing-parameter factor and another factor
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.group.med = "Tuning", 
                                facet.x.med = c("Smarthouse"), 
                                facet.y.med = c("Salinity","Type"),
                                propn.types.med = 0.025))
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.group.med = "Tuning", 
                                facet.x.med = c("Salinity","Type"), 
                                facet.y.med = "Smarthouse",
                                propn.note.med = FALSE))
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.group.med = "Tuning", 
                                facet.x.med = c("Salinity","Type"), 
                                facet.y.med = "Smarthouse",
                                propn.types.med = 0.025))
  testthat::expect_equal(length(t$plots), 1)
  testthat::expect_equal(nrow(t$med.devn.dat), 280)
  testthat::expect_equal(ncol(t$med.devn.dat), 6)
  testthat::expect_true(all(c("Type","Smarthouse", "Salinity", "SmoothParams", "PSA.devn", "DAST") 
                            %in% names(t$med.devn.dat)))
  
    
  #Use plot,by.med
  testthat::expect_silent(
    t <- plotSmoothsMedianDevns(data = smth, response = "PSA", response.smoothed = "sPSA",
                                times = "DAST", trait.types = "response", 
                                plots.by.med = c("Method","Type"), 
                                plots.group.med = c("Tuning"),
                                facet.y.med = c("Smarthouse","Salinity"),
                                propn.types.med = 0.025))
  testthat::expect_equal(length(t$plots), 1)
  testthat::expect_equal(nrow(t$med.devn.dat), 280)
  testthat::expect_equal(ncol(t$med.devn.dat), 6)
  testthat::expect_true(all(c("fac.by", "Smarthouse", "Salinity", "SmoothParams", "PSA.devn", "DAST") 
                            %in% names(t$med.devn.dat)))
  
  testthat::expect_warning(
    smth <- probeSmooths(data = RicePrepped.dat, response = "PSA", response.smoothed = "sPSA", 
                         times = "DAST", 
                         spline.types = c("NCSS", "PS"), smoothing.methods = "logarithmic",
                         df = c(4,7), lambdas = c(0.1,1,10), 
                         which.plots = "medians.dev", 
                         plots.by.med = "Type", 
                         plots.group.med = "Tuning", 
                         facet.y.med = c("Smarthouse","Salinity")),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(nrow(smth), 73920)
  testthat::expect_equal(ncol(smth), 16)
  
  #Test for segmented smoothing
  lambdas <- round(10^c(-0.5, 0, 0.5, 1), digits = 3)
  df <- 4:5
  traits <- c("PSA","PSA.AGR","PSA.RGR")
  DAST.segs <- list(c(-1,5),c(6,13))
  tmp <- subset(RicePrepped.dat, Smarthouse == "NW")
  testthat::expect_warning(
    longiseg.smth <- probeSmooths(data = tmp, 
                                  response = "PSA", response.smoothed = "sPSA",
                                  times = "DAST", 
                                  get.rates = TRUE, 
                                  smoothing.segments = DAST.segs,
                                  spline.types = c("NCSS","PS"), 
                                  df = df, lambdas = list(PS = lambdas), 
                                  npspline.segments = c(4,6), smoothing.methods = "log",
                                  which.plots = "none",
                                  keep.columns = c("Smarthouse", "Salinity")),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(nrow(longiseg.smth), 44352)
  testthat::expect_equal(ncol(longiseg.smth), 16)
  testthat::expect_true(all(unique(longiseg.smth$DAST) == c(-1,1:13)))
  testthat::expect_true(all(c("Type", "TunePar", "TuneVal", "Tuning", "Method", 
                              "Snapshot.ID.Tag", "DAST", "Smarthouse", "Salinity", 
                              "PSA", "PSA.AGR", "PSA.RGR", "sPSA", "sPSA.AGR", "sPSA.RGR") 
                            %in% names(longiseg.smth)))
  #check sPSA.AGR for first DAST of first segment is always NA, 
  #but that first DAST of 2nd segment (DAST 6) is has values that are not NA 
  testthat::expect_true(all(is.na(longiseg.smth$sPSA.AGR[longiseg.smth$DAST == -1]))) 
  testthat::expect_false(all(is.na(longiseg.smth$sPSA.AGR[longiseg.smth$DAST == 6]))) 
  
  #Recalculate using derivatives
  testthat::expect_warning(
    longisegder.smth <- probeSmooths(data = tmp, 
                                     times = "DAST", 
                                     response = "PSA", response.smoothed = "sPSA",
                                     get.rates = TRUE, rates.method = "deriv", 
                                     smoothing.segments = DAST.segs,
                                     spline.types = c("NCSS","PS"), 
                                     df = df, lambdas = list(PS = lambdas), 
                                     npspline.segments = c(4,6), smoothing.methods = "log",
                                     which.plots = "none",
                                     keep.columns = c("Smarthouse", "Salinity")),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(nrow(longisegder.smth), 44352)
  testthat::expect_equal(ncol(longisegder.smth), 16)
  testthat::expect_true(all(unique(longisegder.smth$DAST) == c(-1,1:13)))
  testthat::expect_true(all(c("Type", "TunePar", "TuneVal", "Tuning", "Method", 
                              "Snapshot.ID.Tag", "DAST", "Smarthouse", "Salinity", 
                              "PSA", "PSA.AGR", "PSA.RGR", "sPSA", "sPSA.AGR", "sPSA.RGR") 
                            %in% names(longisegder.smth)))
  #check that, for sPSA.AGR, the DASTs that are missing for differences are not the same as for derivatives
  testthat::expect_true(sum(is.na(longiseg.smth$sPSA.AGR[longiseg.smth$DAST == -1])) !=
                          sum(is.na(longisegder.smth$sPSA[longisegder.smth$DAST == -1]))) 
  #Check that the DAST at the start of each segment has values that are not NA for derivatives 
  testthat::expect_true(!all(is.na(longisegder.smth$sPSA.AGR[longisegder.smth$DAST == -1]))) 
  testthat::expect_true(!all(is.na(longisegder.smth$sPSA.AGR[longisegder.smth$DAST == 6]))) 
  
  #Plot the median deviations
  testthat::expect_warning(
    t <- plotSmoothsMedianDevns(data = longisegder.smth, response = "PSA", 
                                response.smoothed = "sPSA",
                                times = "DAST", trait.types = "AGR", 
                                plots.group.med = "Tuning", 
                                facet.y.med = c("Smarthouse","Salinity"),
                                propn.types.med = 0.025))
  testthat::expect_equal(length(t$plots), 1)
  testthat::expect_equal(nrow(t$med.devn.dat), 168)
  testthat::expect_equal(ncol(t$med.devn.dat), 5)
  #every first DAST must be NA for deviations because PSA.AGR has to be calculated by differences
  testthat::expect_true(all(is.na(t$med.devn.dat$PSA.AGR[t$med.devn.dat$DASTs == -1]))) 
  #every DAST 6, 7 and 13 must be NA for deviations because PSA.AGR has to be calculated 
  #within segments and by differences with ntimes2span == 3 to centre it properly
  testthat::expect_true(all(is.na(t$med.devn.dat$PSA.AGR[t$med.devn.dat$DASTs == 6]))) 
  testthat::expect_true(all(is.na(t$med.devn.dat$PSA.AGR[t$med.devn.dat$DASTs == 7]))) 
  testthat::expect_true(all(is.na(t$med.devn.dat$PSA.AGR[t$med.devn.dat$DASTs == 13]))) 
  testthat::expect_true(all(c("Smarthouse", "Salinity", "SmoothParams", "PSA.AGR.devn", "DAST") 
                            %in% names(t$med.devn.dat)))
  
  #test segmented smooth that is a subset
  lambdas <- round(10^c(-0.5, 0, 0.5, 1), digits = 3)
  df <- 4:5
  traits <- c("PSA","PSA.AGR","PSA.RGR")
  DAST.segs <- list(c(-1,5),c(6,11))
  tmp <- subset(RicePrepped.dat, Smarthouse == "NW")
  testthat::expect_warning(
    longiseg.smth <- probeSmooths(data = tmp, 
                                  response = "PSA", response.smoothed = "sPSA",
                                  times = "DAST", 
                                  get.rates = TRUE, 
                                  smoothing.segments = DAST.segs,
                                  spline.types = c("NCSS","PS"), 
                                  df = df, lambdas = list(PS = lambdas), 
                                  npspline.segments = c(4,6), smoothing.methods = "log",
                                  which.plots = "none",
                                  keep.columns = c("Smarthouse", "Salinity")),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(nrow(longiseg.smth), 38016)
  testthat::expect_equal(ncol(longiseg.smth), 16)
  testthat::expect_true(all(unique(longiseg.smth$DAST) == c(-1,1:11)))
  testthat::expect_true(all(c("Type", "TunePar", "TuneVal", "Tuning", "Method", 
                              "Snapshot.ID.Tag", "DAST", "Smarthouse", "Salinity", 
                              "PSA", "PSA.AGR", "PSA.RGR", "sPSA", "sPSA.AGR", "sPSA.RGR") 
                            %in% names(longiseg.smth)))
  #check sPSA.AGR for first DAST of first segment is always NA, 
  #but that first DAST of 2nd segment (DAST 6) is has values that are not NA 
  testthat::expect_true(all(is.na(longiseg.smth$sPSA.AGR[longiseg.smth$DAST == -1]))) 
  testthat::expect_false(all(is.na(longiseg.smth$sPSA.AGR[longiseg.smth$DAST == 6]))) 
  
  #Recalculate using derivatives
  testthat::expect_warning(
    longisegder.smth <- probeSmooths(data = tmp, 
                                     times = "DAST", 
                                     response = "PSA", response.smoothed = "sPSA",
                                     get.rates = TRUE, rates.method = "deriv", 
                                     smoothing.segments = DAST.segs,
                                     spline.types = c("NCSS","PS"), 
                                     df = df, lambdas = list(PS = lambdas), 
                                     npspline.segments = c(4,6), smoothing.methods = "log",
                                     which.plots = "none",
                                     keep.columns = c("Smarthouse", "Salinity")),
    regexp = "Need at least 4 distinct x values to fit a spline - all fitted values set to NA")
  testthat::expect_equal(nrow(longisegder.smth), 38016)
  testthat::expect_equal(ncol(longisegder.smth), 16)
  testthat::expect_true(all(unique(longisegder.smth$DAST) == c(-1,1:11)))
  testthat::expect_true(all(c("Type", "TunePar", "TuneVal", "Tuning", "Method", 
                              "Snapshot.ID.Tag", "DAST", "Smarthouse", "Salinity", 
                              "PSA", "PSA.AGR", "PSA.RGR", "sPSA", "sPSA.AGR", "sPSA.RGR") 
                            %in% names(longisegder.smth)))
  #check that, for sPSA.AGR, the DASTs that are missing for differences are not the same as for derivatives
  testthat::expect_true(sum(is.na(longiseg.smth$sPSA.AGR[longiseg.smth$DAST == -1])) !=
                          sum(is.na(longisegder.smth$sPSA[longisegder.smth$DAST == -1]))) 
  #every DAST 6, 7 and 13 must be NA for deviations because PSA.AGR has to be calculated 
  #within segments and by differences with ntimes2span == 3 to centre it properly
  testthat::expect_true(!all(is.na(longisegder.smth$sPSA.AGR[longisegder.smth$DAST == -1]))) 
  testthat::expect_true(!all(is.na(longisegder.smth$sPSA.AGR[longisegder.smth$DAST == 6]))) 
  testthat::expect_true(!all(is.na(longisegder.smth$sPSA.AGR[longisegder.smth$DAST == 7]))) 
  testthat::expect_true(!all(is.na(longisegder.smth$sPSA.AGR[longisegder.smth$DAST == 11]))) 
  
})
