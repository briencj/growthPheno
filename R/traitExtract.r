

#Function to extract features from a smoothed trait involving imaging or water use responses
traitExtractFeatures <- function(data, individuals = "Snapshot.ID.Tag", times = "DAP", 
                                 starts.intvl = NULL, stops.intvl = NULL, sep.intvl = "to", 
                                 responses.singletimes = NULL, 
                                 responses.rates = NULL, rates.method = "differences", 
                                 growth.rates = NULL, suffices.rates = NULL, 
                                 water.use = NULL, responses.water = NULL, 
                                 responses.total = NULL, suffix.total = NULL, 
                                 responses.max = NULL, 
                                 times.whole = NULL, 
                                 mergedata = NULL, ...)
{
  #Check options
  options <- c("differences","ratesaverages")
  ratemeth.opt <- options[check.arg.values(rates.method, options=options)]
  options <- c("AGR", "RGR")
  grates <- options[unlist(lapply(growth.rates, check.arg.values, options=options))]
  
  if (is.allnull(times.whole) && !is.allnull(c(starts.intvl, stops.intvl)))
  {
    times.whole <- c(starts.intvl, stops.intvl)
    times.whole <- range(times.whole, na.rm = TRUE)
  }
  
  if (!is.allnull(starts.intvl) && !is.allnull(stops.intvl))
  {
    if (length(starts.intvl) != length(stops.intvl))
      stop("The length of starts.intvl and stops.intvl must be the same")
    if (!all(stops.intvl > starts.intvl))
      stop("Each element of stops.intvls should be the greater than the corresponding element of starts.intvl")
  }
  
  if (!is.allnull(growth.rates) && !is.allnull(suffices.rates) && length(growth.rates) != suffices.rates)
    stop("The length of growth.rates and suffices.rates should be the same")
  
  #Get a complete set of IDs 
  indv.dat <- data.frame(sort(unique(data[[individuals]])))
  names(indv.dat) <- individuals
  
  ####Get the single times
  if (!is.allnull(responses.singletimes))
  {
    times1 <- sort(unique(c(starts.intvl, stops.intvl)))
    if (is.allnull(times1))
      stop("No times have been specified using starts.intvl or stops.intvl")
    for (t1 in times1)
      indv.dat <- merge(indv.dat, 
                        getTimesSubset(data = data, 
                                       responses = responses.singletimes, 
                                       individuals = individuals, times = times, 
                                       which.times = t1, 
                                       suffix = t1, 
                                       include.individuals = TRUE),
                        by = individuals, sort = FALSE)
  }
  
  ####Get the growth rates
  if (!is.allnull(responses.rates))
  {
    if (is.allnull(suffices.rates))
      suffice.rates <- grates
    suffices.intvl <- paste(starts.intvl, stops.intvl, sep = sep.intvl)
    
    if (ratemeth.opt == "differences")
    {  ### Rates for specific intervals from the smoothed data by differencing
      for (r in responses.rates)
      { 
        for (k in 1:length(suffices.intvl))
        { 
          indv.dat <- merge(indv.dat, 
                            byIndv4Intvl_GRsDiff(data = data, responses = r, 
                                                 individuals = individuals, times = times, 
                                                 which.rates = grates,
                                                 suffices.rates = suffice.rates, 
                                                 start.time = starts.intvl[k], 
                                                 end.time = stops.intvl[k], 
                                                 suffix.interval = suffices.intvl[k]),
                            by = individuals, sort = FALSE)
        }
      }
    } else # derivatives
    {
      if (!all(unlist(lapply(grates, function(x, responses) paste(responses, x, sep = "."), responses = responses.rates)) 
               %in% names(data)))
        stop("There is not a column in data for every combinations of reponses.rates and growth rates")
      for (r in responses.rates)
      { 
        for (k in 1:length(suffices.intvl))
        { 
          indv.dat <- merge(indv.dat, 
                            byIndv4Intvl_GRsAvg(data = data, responses = r, 
                                                individuals = individuals, times = times, 
                                                which.rates = grates,
                                                suffices.rates = suffice.rates, 
                                                start.time = starts.intvl[k], 
                                                end.time = stops.intvl[k], 
                                                suffix.interval = suffices.intvl[k],
                                                sep = "_"),
                            by = individuals, sort = FALSE)
        }
      }
    }
  }
  
  ### Get the water traits
  if (!is.null(water.use) && !is.allnull(responses.water))
  {
    suffix.rate <- "R"
    if (grepl(".", water.use, fixed = TRUE))
      suffix.rate <- ".Rate"
    suffix.index <- "I"
    if (grepl(".", water.use, fixed = TRUE))
      suffix.index <- ".Index"
    for (r in responses.rates)
    { 
      for (k in 1:length(suffices.intvl))
      { 
        indv.dat <- merge(indv.dat, 
                      byIndv4Intvl_WaterUse(data = data, 
                                            water.use = water.use, responses = r, 
                                            individuals = individuals, times = times, 
                                            trait.types =c("WU", "WUR", "WUI"),
                                            suffix.rate = suffix.rate, 
                                            suffix.index = suffix.index, 
                                            start.time = starts.intvl[k], 
                                            end.time = stops.intvl[k], 
                                            suffix.interval = suffices.intvl[k]),
                      by = individuals, sort = FALSE)
      }
    }
  }
 
  #get the total over the whole imaging period
  if (!is.allnull(responses.total))
  {
    if (is.null(suffix.total))
      suffix.total = paste(times.whole[1], times.whole[2], sep = sep.intvl)
    for (r in responses.total)
    { 
      indv.dat <- merge(indv.dat, 
                        byIndv4Intvl_ValueCalc(data = data, 
                                               response = r, 
                                               individuals = individuals, times = times, 
                                               FUN = "sum", addFUN2name = FALSE, 
                                               start.time = times.whole[1], 
                                               end.time = times.whole[2], 
                                               suffix.interval = suffix.total))
    }
  }
  
  #get max over the whole imaging period and the time at which it occurred
  if (!is.allnull(responses.max))
  {
    for (r in responses.max)
    { 
      indv.dat <- merge(indv.dat, 
                        byIndv4Intvl_ValueCalc(data = data, 
                                               response = r, 
                                               individuals = individuals, times = times, 
                                               FUN = "max", which.obs = FALSE, which.values = times, 
                                               start.time = times.whole[1], 
                                               end.time = times.whole[2], 
                                               suffix.interval = NULL),
                        by = individuals, sort = FALSE)
    }
  }

  if (!is.null(mergedata))
    indv.dat <- merge(mergedata, indv.dat, sort = FALSE)
    
  return(indv.dat)
}