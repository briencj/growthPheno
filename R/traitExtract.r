

#Function to extract features from a smoothed trait involving imaging or water use responses
traitExtractFeatures <- function(data, individuals = "Snapshot.ID.Tag", times = "DAP", 
                                 starts.intvl = NULL, stops.intvl = NULL, 
                                 suffices.intvl = NULL, 
                                 responses4intvl.rates = NULL, 
                                 growth.rates = NULL, 
                                 growth.rates.method = "differences", 
                                 suffices.growth.rates = NULL, 
                                 water.use4intvl.traits = NULL, responses4water = NULL, 
                                 water.trait.types = c("WU", "WUR", "WUI"), 
                                 suffix.water.rate = "R", suffix.water.index = "I", 
                                 responses4singletimes = NULL, times.single = NULL, 
                                 responses4overall.rates = NULL, 
                                 water.use4overall.water = NULL, responses4overall.water = NULL, 
                                 responses4overall.totals = NULL, 
                                 responses4overall.max = NULL, 
                                 intvl.overall = NULL, suffix.overall = NULL, 
                                 sep.times.intvl = "to", sep.suffix.times = ".", 
                                 sep.growth.rates = ".", sep.water.traits = "", 
                                 mergedata = NULL, ...)
{
  #Check that specified columns are in data
  req.cols <- unique(c(individuals, times, responses4singletimes, 
                       responses4intvl.rates, responses4overall.rates,
                       water.use4intvl.traits, responses4water, water.use4overall.water, responses4overall.water,
                       responses4overall.totals, responses4overall.max))
  checkNamesInData(req.cols, data)
    
  #Check options
  options <- c("differences","ratesaverages")
  respratemeth.opt <- options[check.arg.values(growth.rates.method, options=options)]
  options <- c("AGR", "RGR")
  grates <- options[unlist(lapply(growth.rates, check.arg.values, options=options))]
  
  options <- c("WU", "WUR", "WUI")
  water.traits <- options[unlist(lapply(water.trait.types, check.arg.values, options=options))]
  water.traits <- c("WU", "WUR", "WUI")[c("WU", "WUR", "WUI") %in% water.traits]
  names(water.traits) <- water.traits
  
  #Process times arguments
  times.all <- NULL
  if(!is.allnull(c(starts.intvl, stops.intvl)))
    times.all <- c(starts.intvl, stops.intvl)
  #Process intvl.overall
  if (!is.allnull(intvl.overall))
  {
    if (length(intvl.overall) !=  2 && intvl.overall[1] >= intvl.overall[1])
      stop("intvl.overall should contain the two values, the second value being greater than the first")
  }
  else 
  {
    if(!is.null(times.all))
      intvl.overall <- range(times.all, na.rm = TRUE)
  }
  
  #Process times.single
  if (is.null(times.single)  && !is.null(times.all))
    times.single <- sort(unique(c(starts.intvl, stops.intvl)))
  
  #Check starts.intvl and stops.intvl and set up suffices.intvl
  if (is.allnull(starts.intvl) || is.allnull(stops.intvl))
  {   
    if (!is.allnull(c(responses4intvl.rates, water.use4intvl.traits, responses4water)))
      stop("Neither starts.intvl nor stops.intvl can be NULL")
  }
  else
  {
    if (length(starts.intvl) != length(stops.intvl))
      stop("The length of starts.intvl and stops.intvl must be the same")
    if (!all(stops.intvl > starts.intvl))
      stop("Each element of stops.intvls should be the greater than the corresponding element of starts.intvl")
    if (is.null(suffices.intvl))
    { 
      if (is.null(sep.times.intvl))
      {
        if (length(starts.intvl) != 1)
          stops.intvl(paste0("Both suffices.intvl and sep.times.intvl are NULL, ",
                             "start.intvl and stop.intvl specify more than one interval "))
      } else
        suffices.intvl <- paste(starts.intvl, stops.intvl, sep = sep.times.intvl)
    } else #use suffices.intvl
    {
      if (length(suffices.intvl) != length(starts.intvl))
        stop("The length of suffices.intvl does not match starts.intvl and stops.intvl")
    }
  }
  
  if (!is.allnull(growth.rates) && !is.allnull(suffices.growth.rates) && length(growth.rates) != suffices.growth.rates)
    stop("The length of growth.rates and suffices.growth.rates should be the same")
  
  #Get a complete set of IDs 
  indv.dat <- data.frame(sort(unique(data[[individuals]])))
  names(indv.dat) <- individuals
  
  ####Get the single times
  if (!is.allnull(responses4singletimes))
  {
    if (is.allnull(times.single))
      stop("No times available for single-time traits")
    for (t1 in times.single)
      indv.dat <- merge(indv.dat, 
                        getTimesSubset(data = data, 
                                       responses = responses4singletimes, 
                                       individuals = individuals, times = times, 
                                       which.times = t1, 
                                       suffix = t1, 
                                       sep.suffix.times = sep.suffix.times,
                                       include.individuals = TRUE),
                        by = individuals, sort = FALSE)
  }
  
  ####Get the growth rates
  if (!is.allnull(responses4intvl.rates))
  {
    if (length(grates) == 0)
      stop("growth.rates needs to be set for responses4intvl.rates")
    if (is.allnull(suffices.growth.rates))
      suffices.growth.rates <- grates
    
    if (respratemeth.opt == "differences")
    {  ### Rates for specific intervals from the smoothed data by differencing
      for (r in responses4intvl.rates)
      { 
        for (k in 1:length(suffices.intvl))
        { 
          indv.dat <- merge(indv.dat, 
                            byIndv4Intvl_GRsDiff(data = data, responses = r, 
                                                 individuals = individuals, times = times, 
                                                 which.rates = grates,
                                                 suffices.rates = suffices.growth.rates, 
                                                 sep.rates = sep.growth.rates, 
                                                 start.time = starts.intvl[k], 
                                                 end.time = stops.intvl[k], 
                                                 suffix.interval = suffices.intvl[k], 
                                                 sep.suffix.interval = sep.suffix.times),
                            by = individuals, sort = FALSE)
        }
      }
    } else # derivatives
    {
      req.cols <- unlist(lapply(grates, 
                                function(x, responses) paste(responses, x, sep = sep.growth.rates), 
                                responses = responses4intvl.rates))
      if (!all(req.cols %in% names(data)))
        stop(paste(req.cols[!(req.cols %in% names(data))], collapse = ", "), 
             " is/are not in data as is required for the ratesaverages option")
      for (r in responses4intvl.rates)
      { 
        for (k in 1:length(suffices.intvl))
        { 
          indv.dat <- merge(indv.dat, 
                            byIndv4Intvl_GRsAvg(data = data, responses = r, 
                                                individuals = individuals, times = times, 
                                                which.rates = grates,
                                                suffices.rates = suffices.growth.rates, 
                                                sep.rates = sep.growth.rates, 
                                                start.time = starts.intvl[k], 
                                                end.time = stops.intvl[k], 
                                                suffix.interval = suffices.intvl[k], 
                                                sep.suffix.interval = sep.suffix.times,
                                                sep.levels = "_"),
                            by = individuals, sort = FALSE)
        }
      }
    }
  }
  
  ### Get the water traits
  if (!is.allnull(water.use4intvl.traits) && !is.allnull(responses4water))
  {
    # suffix.rate <- "R"
    # if (grepl(".", water.use4intvl.traits, fixed = TRUE))
    #   suffix.rate <- ".Rate"
    # suffix.index <- "I"
    # if (grepl(".", water.use4intvl.traits, fixed = TRUE))
    #   suffix.index <- ".Index"
    if (length(water.use4intvl.traits) > 1 && length(responses4water) > 1)
    {
      if (length(water.use4intvl.traits) != length(responses4water))
        stop(paste0("If both water.use4intvl.traits and responses4water contain multiple values, the number must be the same so that ",
                    "they can be processed in parallel"))
      #Process water.use4intvl.traits and responses4water in parallel
      for (i in 1:length(water.use4intvl.traits))
      { 
        for (k in 1:length(suffices.intvl))
        { 
          indv.dat <- merge(indv.dat, 
                            byIndv4Intvl_WaterUse(data = data, 
                                                  water.use = water.use4intvl.traits[i], 
                                                  responses = responses4water[i], 
                                                  individuals = individuals, times = times, 
                                                  trait.types = water.traits,
                                                  suffix.rate = suffix.water.rate, 
                                                  suffix.index = suffix.water.index, 
                                                  sep.water.traits = sep.water.traits, 
                                                  sep.responses = sep.growth.rates, 
                                                  start.time = starts.intvl[k], 
                                                  end.time = stops.intvl[k], 
                                                  suffix.interval = suffices.intvl[k], 
                                                  sep.suffix.interval = sep.suffix.times),
                            by = individuals, sort = FALSE)
        }
      }
    } else
    {
      #Process water.use4intvl.traits and responses4water in combination
      for (w in water.use4intvl.traits)
      {     
        for (r in responses4water)
        { 
          for (k in 1:length(suffices.intvl))
          { 
            indv.dat <- merge(indv.dat, 
                              byIndv4Intvl_WaterUse(data = data, 
                                                    water.use = w, responses = r, 
                                                    individuals = individuals, times = times, 
                                                    trait.types = water.traits,
                                                    suffix.rate = suffix.water.rate, 
                                                    suffix.index = suffix.water.index, 
                                                    sep.water.traits = sep.water.traits, 
                                                    sep.responses = sep.growth.rates, 
                                                    start.time = starts.intvl[k], 
                                                    end.time = stops.intvl[k], 
                                                    suffix.interval = suffices.intvl[k], 
                                                    sep.suffix.interval = sep.suffix.times),
                              by = individuals, sort = FALSE)
          }
        }
      }
    }
  }
  
  #Deal with the overall traits
  ####Get the overall responses rates
  if (!is.allnull(responses4overall.rates))
  {
    if (is.allnull(intvl.overall))
      stop("No times available for overall traits")
    if (length(grates) == 0)
      stop("growth.rates needs to be set for responses4overall.rates")
    if (is.allnull(suffices.growth.rates))
      suffices.growth.rates <- grates

    if (respratemeth.opt == "differences")
    {  ### Rates for specific intervals from the smoothed data by differencing
      for (r in responses4overall.rates)
      { 
        indv.dat <- merge(indv.dat, 
                          byIndv4Intvl_GRsDiff(data = data, responses = r, 
                                               individuals = individuals, times = times, 
                                               which.rates = grates, 
                                               suffices.rates = suffices.growth.rates, 
                                               sep.rates = sep.growth.rates, 
                                               start.time = intvl.overall[1], 
                                               end.time = intvl.overall[2], 
                                               suffix.interval = suffix.overall, 
                                               sep.suffix.interval = sep.suffix.times),
                          by = individuals, sort = FALSE)
      }
    } else # derivatives
    {
      req.cols <- unlist(lapply(grates, 
                                function(x, responses) paste(responses, x, sep = sep.growth.rates), 
                                responses = responses4overall.rates))
      if (!all(req.cols %in% names(data)))
        stop(paste(req.cols[!(req.cols %in% names(data))], collapse = ", "), 
             " is/are not in data as is required for the ratesaverages option")
      for (r in responses4overall.rates)
      { 
        indv.dat <- merge(indv.dat, 
                          byIndv4Intvl_GRsAvg(data = data, responses = r, 
                                              individuals = individuals, times = times, 
                                              which.rates = grates,
                                              suffices.rates = suffices.growth.rates, 
                                              sep.rates = sep.growth.rates, 
                                              start.time = intvl.overall[1], 
                                              end.time = intvl.overall[2], 
                                              suffix.interval = suffix.overall, 
                                              sep.suffix.interval = sep.suffix.times,
                                              sep.levels = "_"),
                          by = individuals, sort = FALSE)
      }
    }
  }
  
  
  ### Get the overall water traits
  if (!is.allnull(water.use4overall.water) && !is.allnull(responses4overall.water))
  {
    # suffix.rate <- "R"
    # if (grepl(".", water.use, fixed = TRUE))
    #   suffix.rate <- ".Rate"
    # suffix.index <- "I"
    # if (grepl(".", water.use, fixed = TRUE))
    #   suffix.index <- ".Index"
    if (is.allnull(intvl.overall))
      stop("No times available for overall traits")
    if (length(water.use4overall.water) > 1 && length(responses4overall.water) > 1)
    {
      if (length(water.use4overall.water) != length(responses4overall.water))
        stop(paste0("If both water.use4overall.water and overallresponses.water contain multiple values, the number must be the same so that ",
                    "they can be processed in parallel"))
      #Process water.use and responses.water in parallel
      for (i in 1:length(water.use4overall.water))
      { 
        indv.dat <- merge(indv.dat, 
                          byIndv4Intvl_WaterUse(data = data, 
                                                water.use = water.use4overall.water[i], responses = responses4overall.water[i], 
                                                individuals = individuals, times = times, 
                                                trait.types = water.traits,
                                                suffix.rate = suffix.water.rate, 
                                                suffix.index = suffix.water.index, 
                                                sep.water.traits = sep.water.traits, 
                                                sep.responses = sep.growth.rates, 
                                                start.time = intvl.overall[1], 
                                                end.time = intvl.overall[2], 
                                                suffix.interval = suffix.overall, 
                                                sep.suffix.interval = sep.suffix.times),
                          by = individuals, sort = FALSE)
      }
    } else
    {
      #Process water.use4overall.water and responses4overall.water in combination
      for (w in water.use4overall.water)
      {     
        for (r in responses4overall.water)
        { 
          indv.dat <- merge(indv.dat, 
                            byIndv4Intvl_WaterUse(data = data, 
                                                  water.use = w, responses = r, 
                                                  individuals = individuals, times = times, 
                                                  trait.types = water.traits,
                                                  suffix.rate = suffix.water.rate, 
                                                  sep.water.traits = sep.water.traits, 
                                                  suffix.index = suffix.water.index, 
                                                  start.time = intvl.overall[1], 
                                                  end.time = intvl.overall[2], 
                                                  suffix.interval = suffix.overall, 
                                                  sep.suffix.interval = sep.suffix.times),
                            by = individuals, sort = FALSE)
        }
      }
    }
  }
  
  #get totals over the whole imaging period
  if (!is.allnull(responses4overall.totals))
  {
    for (r in responses4overall.totals)
    { 
      indv.dat <- merge(indv.dat, 
                        byIndv4Intvl_ValueCalc(data = data, 
                                               response = r, 
                                               individuals = individuals, times = times, 
                                               FUN = "sum", addFUN2name = FALSE, 
                                               start.time = intvl.overall[1], 
                                               end.time = intvl.overall[2], 
                                               suffix.interval = suffix.overall, 
                                               sep.suffix.interval = sep.suffix.times))
    }
  }
  
  #get max over the whole imaging period and the time at which it occurred
  if (!is.allnull(responses4overall.max))
  {
    for (r in responses4overall.max)
    { 
      indv.dat <- merge(indv.dat, 
                        byIndv4Intvl_ValueCalc(data = data, 
                                               response = r, 
                                               individuals = individuals, times = times, 
                                               FUN = "max", which.obs = FALSE, which.values = times, 
                                               start.time = intvl.overall[1], 
                                               end.time = intvl.overall[2], 
                                               suffix.interval = NULL),
                        by = individuals, sort = FALSE)
    }
  }

  if (!is.null(mergedata))
    indv.dat <- merge(mergedata, indv.dat, sort = FALSE)
    
  return(indv.dat)
}