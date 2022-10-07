#Functions to produce a value in an interval for each interval
#Function to form the growth rates over an interval for a set of responses
"byIndv4Intvl_GRsDiff" <- function(data, responses, 
                                   individuals = "Snapshot.ID.Tag", times = "DAP", 
                                   which.rates = c("AGR","PGR","RGR"), 
                                   suffices.rates=NULL, sep.rates = ".", 
                                   start.time, end.time, 
                                   suffix.interval, sep.suffix.interval = ".")
{ 
  options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) && length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  
  if (!all(sapply(list(start.time, end.time), 
                  function(x) length(x) == 1)))
    stop("At least one of start.time and end.time are not length one")
  if (!is.allnull(suffix.interval))
  { 
    if (length(suffix.interval) == 1)
      suffix.interval <- paste0(sep.suffix.interval, suffix.interval) #so can use paste0, which allows for NULL suffix.interval
    else
      stop("suffix.interval is not length one")
  }
  
  if (is.null(sep.rates))
    stop('sep.rates cannot be NULL; for no separator use "."')
  
  #Check that responses, individuals and times are in data
  checkNamesInData(c(responses, individuals, times), data = data)
  
  interval.resp <- mapply(function(time, suffix, data, responses, times) 
                          {
                            dat <- getTimesSubset(data = data, responses = responses, 
                                           times = times, which.times = time, 
                                           include.times = TRUE, 
                                           suffix = suffix)
                            names(dat)[grepl(individuals, names(dat))] <- individuals
                            return(dat)
                          }, time = c(start.time,end.time), suffix = c("start", "end"), 
                          MoreArgs = list(data = data, responses = c(individuals, responses), 
                                             times = times), SIMPLIFY = FALSE)
  interval.resp <- merge(interval.resp[[1]], interval.resp[[2]], all = TRUE, sort = FALSE)
                          
  times.start <- paste(times,"start",sep=".")
  times.end <- paste(times,"end",sep=".")
  
  interval.resp[times.start] <- convertTimes2numeric(interval.resp[[times.start]])
  interval.resp[times.end] <- convertTimes2numeric(interval.resp[[times.end]])
  growth.rates <- lapply(responses, 
                         function(name, interval.resp, start.time, end.time, 
                                  which.rates, times = "DAP")
                         { #check which.rates
                           if (any(is.na(pmatch(which.rates, c("AGR","PGR","RGR")))))
                             stop("which.rates has at least one illegal value")
                           start.name <- paste(name,"start",sep=".")
                           start.time <- paste(times,"start",sep=".")
                           end.name <- paste(name,"end",sep=".")
                           end.time <- paste(times,"end",sep=".")
                           rates <- vector("list", length = length(which.rates))
                           k <- 0
                           if ("AGR" %in% which.rates)
                           { 
                             k <- k + 1
                             AGR <- (interval.resp[end.name] - interval.resp[start.name]) / 
                               (interval.resp[end.time] - interval.resp[start.time])
                             rates[k] <- AGR
                             if (is.null(suffices.rates))
                               names(rates)[k] <- paste0(name, sep.rates, "AGR", suffix.interval)
                             else
                               names(rates)[k] <- paste0(name, sep.rates, suffices.rates[match("AGR",opt)],
                                                         suffix.interval)
                           }
                           if ("RGR" %in% which.rates)
                           { 
                             k <- k + 1
                             RGR <- (log(interval.resp[end.name]) - log(interval.resp[start.name])) / 
                               (interval.resp[end.time] - interval.resp[start.time])
                             rates[k] <- RGR
                             if (is.null(suffices.rates))
                               names(rates)[k] <- paste0(name, sep.rates, "RGR", suffix.interval)
                             else
                               names(rates)[k] <- paste0(name, sep.rates, suffices.rates[match("RGR",opt)],
                                                         suffix.interval)
                           }
                           if ("PGR" %in% which.rates)
                           { 
                             k <- k +1
                             if ("RGR" %in% which.rates)
                               PGR <- exp(RGR)
                             else
                             { 
                               PGR <- (log(interval.resp[end.name]) - log(interval.resp[start.name])) / 
                                 (interval.resp[end.time] - interval.resp[start.time])
                               PGR <- exp(PGR)
                             }
                             rates[k] <- PGR
                             if (is.null(suffices.rates))
                               names(rates)[k] <- paste0(name, sep.rates, "PGR", suffix.interval)
                             else
                               names(rates)[k] <- paste0(name, sep.rates, suffices.rates[match("PGR",opt)],
                                                         suffix.interval)
                           }
                           rates <- as.data.frame(rates)
                           return(rates)
                         },
                         interval.resp = interval.resp, 
                         start.time = start.time, end.time = end.time, 
                         times=times, which.rates = opt)
  growth.rates <- as.data.frame(growth.rates)
  growth.rates <- cbind(interval.resp[individuals], growth.rates)
  rownames(growth.rates) <- NULL
  #use getTimesSubset to add the individuals to the growth.rates data.frame  
  # growth.rates <- merge(growth.rates, 
  #                       getTimesSubset(data = data, responses = individuals, times = times, 
  #                                      which.times = start.time, suffix = NULL))
  return(growth.rates)
}

#Function to form the growth rates over an interval for a set of responses by averaging growth rates
"byIndv4Intvl_GRsAvg" <- function(data, responses, 
                                  individuals = "Snapshot.ID.Tag", times = "DAP", 
                                  which.rates = c("AGR","RGR"), 
                                  suffices.rates=c("AGR","RGR"), sep.rates = ".", 
                                  start.time, end.time, 
                                  suffix.interval, sep.suffix.interval = ".", 
                                  sep.levels=".", na.rm=FALSE)
{  
  #Check that responses, individuals and times are in data
  checkNamesInData(c(responses, individuals, times), data = data)

  if (!all(sapply(list(start.time, end.time), 
           function(x) length(x) == 1)))
    stop("At least one of start.time and end.time are not length one")
  if (!is.allnull(suffix.interval))
  { 
    if (length(suffix.interval) == 1)
      suffix.interval <- paste0(sep.suffix.interval, suffix.interval) #so can use paste0, which allows for NULL suffix.interval
    else
      stop("suffix.interval is not length one")
  }
  
  if (is.null(sep.rates))
    stop('sep.rates cannot be NULL; for no separator use "."')
  if (!(is.allnull(suffices.rates)) && length(which.rates) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  #Check that required growth rates are in data
  if (is.allnull(suffices.rates))
    response.grates <- responses
  else
    response.grates <- unlist(lapply(suffices.rates,
                                   function(grate, responses)
                                   { 
                                     response.grates <- paste(responses, grate, sep = sep.rates)
                                     if (!all(response.grates %in% names(data)))
                                       stop("The following growth rates are not in data: ",
                                            paste0(response.grates[!(response.grates %in% names(data))]))
                                     return(response.grates)
                                   },
                                   responses=responses))
  #Get data for the times 
  times.vals <- unique(convertTimes2numeric(data[[times]]))
  times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
  combos <- expand.grid(sort(unique(data[[individuals]])), times.vals)
  names(combos) <- c(individuals, times)
  interval.resp <- getTimesSubset(data = data, responses = response.grates, 
                                  times = times, which.times = times.vals, 
                                  include.times = TRUE, include.individuals = TRUE)
  interval.resp[times] <- convertTimes2numeric(interval.resp[[times]])
  interval.resp <- merge(combos, interval.resp, all = TRUE, sort = FALSE)
  #calculate the weights
  interval.resp <- split(interval.resp, as.list(interval.resp[individuals]), sep=sep.levels)
  interval.resp <- lapply(interval.resp, 
                          function(data, times = "DAP")
                          { 
                            n <- nrow(data)
                            n1 <- n - 1
                            data[2:n, times] <- data[2:n, times] - 
                              data[1:n1, times]
                            data[1, times] <- data[2, times]
                            data[2:n1, times] <- (data[1:(n1-1), times] + 
                                                    data[2:n1, times])/2
                            return(data)
                          }, times = times)
  interval.resp <- do.call(rbind, interval.resp)
  
  #Calculate growth rates by averaging within an interval
  avgrowth.rates <- lapply(responses, 
                           function(name, interval.resp, which.rates = c("AGR","RGR"))
                           { #check which.rates
#                             if (any(is.na(pmatch(which.rates, c("AGR","RGR")))))
#                               stop("which.rates has at least one illegal value")
                             for (k in 1:length(which.rates))
                             {
                               if (is.allnull(suffices.rates))
                                 name.rate <- name
                               else
                                 name.rate <- paste(name, suffices.rates[k], sep = sep.rates)
                               if (which.rates[k] == "RGR")
                                 interval.resp[name.rate] <- log(interval.resp[name.rate])
                               if (k==1)
                                 rates <- byIndv_ValueCalc(data=interval.resp, response = name.rate, 
                                                           individuals = individuals, 
                                                           FUN = "weighted.mean", 
                                                           weights=times, 
                                                           na.rm=na.rm, sep.levels=sep.levels)
                               else
                                 rates <- merge(rates, 
                                                byIndv_ValueCalc(data=interval.resp, response = name.rate, 
                                                                 individuals = individuals, 
                                                                 FUN = "weighted.mean", 
                                                                 weights=times, 
                                                                 na.rm=na.rm, sep.levels=sep.levels),
                                 )
                               if (which.rates[k] == "RGR")
                                 rates[length(rates)] <- exp(rates[length(rates)])
                               names(rates)[length(rates)] <- paste0(name.rate, suffix.interval)
                             }
                             return(rates)
                           },
                           interval.resp = interval.resp, 
                           which.rates = which.rates)
  avgrowth.rates <- as.data.frame(avgrowth.rates)
  return(avgrowth.rates)
}

#Function to calculate a value for observations within an interval for a set of responses
"byIndv4Intvl_ValueCalc" <- function(data, response, 
                                     individuals = "Snapshot.ID.Tag", times = "DAP", 
                                     FUN = "max", which.obs = FALSE, which.values = NULL, 
                                     addFUN2name = TRUE, sep.FUNname = ".", 
                                     start.time=NULL, end.time=NULL, 
                                     suffix.interval=NULL, sep.suffix.interval = ".",
                                     sep.levels=".", weights=NULL, na.rm=TRUE, ...)
{  
  #Trap which.levels and give message to replace it with which.values is not set
  tempcall <- list(...)
  if (length(tempcall) && "which.levels" %in% names(tempcall))
    stop("replace which.levels with which.values")
  
  #Check that response, individuals and times are in data
  checkNamesInData(c(response, individuals, times), data = data)
  
  if (!all(sapply(list(start.time, end.time, suffix.interval), 
                  function(x) length(x) <= 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length one or NULL")
  
  #Get data for the times
  if (all(is.null(c(start.time, end.time))))
    interval.resp <- data[c(individuals,response,times)]
  else
  { 
    times.vals <- unique(convertTimes2numeric(data[[times]]))
    if (is.null(start.time))
      times.vals <- times.vals[times.vals <= end.time]
    else
    { 
      if (is.null(end.time))
        times.vals <- times.vals[times.vals >= start.time]
      else
        times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
    }
    interval.resp <- getTimesSubset(data = data, responses = response, 
                                    times = times, which.times = times.vals, 
                                    include.times = TRUE)
    interval.resp <- cbind(getTimesSubset(data = data, responses = individuals, 
                                          times = times, which.times = times.vals, 
                                          include.times = FALSE),
                           interval.resp)
  }
  
  #Calculate a value within an interval for each individual
  val.dat <- byIndv_ValueCalc(data = interval.resp, response=response, individuals = individuals, 
                              FUN = FUN, which.obs = which.obs, which.values = which.values, 
                              addFUN2name = addFUN2name, weights=weights, na.rm=na.rm, sep.levels=sep.levels, ...)
  if (!is.null(suffix.interval))
  {
    if (addFUN2name)
    {   
      names(val.dat)[match(paste(response, FUN, sep=sep.FUNname), names(val.dat))] <- 
        paste(paste(response, FUN, sep = sep.FUNname), suffix.interval, sep = sep.suffix.interval)
      if (which.obs)
        names(val.dat)[match(paste(paste(response, FUN, sep = sep.FUNname), "obs", sep="."), names(val.dat))] <- 
          paste(paste(response, FUN, sep = sep.FUNname), "obs", suffix.interval, sep = sep.suffix.interval)
      if (!is.null(which.values))
        names(val.dat)[match(paste(paste(response, FUN, sep = sep.FUNname), which.values, sep="."), names(val.dat))] <- 
          paste(paste(response, FUN, sep = sep.FUNname), which.values, suffix.interval, sep = sep.suffix.interval)
    } else
    {
      names(val.dat)[match(response, names(val.dat))] <- 
        paste(response, suffix.interval, sep = sep.suffix.interval)
      if (which.obs)
        names(val.dat)[match(paste(response, "obs", sep="."), names(val.dat))] <- 
          paste(response, "obs", suffix.interval, sep = sep.suffix.interval)
      if (!is.null(which.values))
        names(val.dat)[match(paste(response, which.values, sep="."), names(val.dat))] <- 
          paste(response, which.values, suffix.interval, sep = sep.suffix.interval)
    }
  }
  
  return(val.dat)
}

# Function to calculate water use indices (WUI) over an interval for a set of responses
"byIndv4Intvl_WaterUse" <- function(data, water.use = "Water.Use", responses = NULL, 
                                    individuals = "Snapshot.ID.Tag", times = "DAP", 
                                    trait.types = c("WU", "WUR", "WUI"), 
                                    suffix.rate = "R", suffix.index = "I", 
                                    sep.water.traits = "", sep.responses = ".", 
                                    start.time, end.time, 
                                    suffix.interval = NULL, sep.suffix.interval = ".", 
                                    na.rm = FALSE)
{ 
  options <- c("WU", "WUR", "AGR", "WUI", "all")
  traits <- options[unlist(lapply(trait.types, check.arg.values, options=options))]
  if ("all" %in% traits)
    traits <- c("WU", "WUR", "WUI")
  traits <- c("WU", "WUR", "WUI")[c("WU", "WUR", "WUI") %in% traits]
  names(traits) <- traits

  if (is.null(responses) && "WUI" %in% traits)
    stop("Must specify a response to calculate a WUI")
  
  #Check that responses, individuals and times are in data
  checkNamesInData(c(responses, water.use, individuals, times), data = data)
  
  if (!all(sapply(list(suffix.rate, suffix.index, start.time, end.time, suffix.interval), 
                  function(x) length(x) <= 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length zero or one")
  
  #get the water use always because needed for WUI
  WU.dat <- byIndv_ValueCalc(data = subset(data, 
                                           convertTimes2numeric(eval(parse(text=times))) >= 
                                             max(start.time)+1 & 
                                             convertTimes2numeric(eval(parse(text=times))) <= 
                                             max(end.time)),
                             response = water.use, individuals = individuals, 
                             FUN = "sum", na.rm = na.rm)
  names(WU.dat)[match(paste(water.use,"sum",sep="."), names(WU.dat))] <- "WU"
  
  #Construct trait.names
  trait.names <- c("WU", "WUR", "WUI")
  names(trait.names) <- trait.names
  trait.names["WU"] <- water.use
  trait.names["WUR"] <- paste(water.use, suffix.rate, sep = sep.water.traits)
  trait.names["WUI"] <- paste(water.use, suffix.index, sep = sep.water.traits)
  trait.names <- trait.names[traits]
  if (length(suffix.interval) == 1)
    trait.names <- paste(trait.names, suffix.interval, 
                         sep = ifelse(is.null(suffix.interval), "", sep.suffix.interval))
  names(trait.names) <- traits

  #Get the values of the responses for the start and end of the time interval
  interval.resp <- cbind(getTimesSubset(data = data, responses = individuals, 
                                        times = times, which.times = start.time, 
                                        suffix = NULL))
  if (!is.null(responses))
    interval.resp <- cbind(interval.resp,
                           getTimesSubset(data = data, responses = responses, 
                                          times = times, which.times = start.time, 
                                          include.times = TRUE, 
                                          suffix = "start"), 
                           getTimesSubset(data = data, responses = responses, 
                                          times = times, which.times = end.time, 
                                          include.times = TRUE, 
                                          suffix = "end"))
  interval.resp <- merge(interval.resp, WU.dat, by = individuals)
  interval.resp <- interval.resp[ do.call(order, interval.resp), ]
  if ("WUR" %in% traits)
    interval.resp["WUR"] <- interval.resp["WU"][,1]/ (end.time - start.time)
  #Reduce WU.dat to the requested WU.traits
  WU.traits <- c("WU","WUR")[c("WU","WUR") %in% traits]
  WU.dat <- interval.resp[c(individuals, WU.traits)]
  names(WU.dat) <- c(individuals, trait.names[WU.traits])
  
  #Get the WUI trait, as requested
  if ("WUI" %in% traits)
  {  
    resp.water <- lapply(responses, 
                         function(name, interval.resp, start.time, end.time, trait.names, 
                                  traits)
                         { 
                           start.name <- paste(name,"start",sep=".")
                           end.name <- paste(name,"end",sep=".")
                           #Calculate the water traits
                           rates <- WUI(interval.resp[end.name][,1] - interval.resp[start.name][,1], 
                                        interval.resp["WU"][,1])
                           return(rates)
                         },
                         interval.resp = interval.resp, 
                         start.time = start.time, end.time = end.time, 
                         trait.names = trait.names, traits = traits)
    names(resp.water) <- paste(responses, trait.names["WUI"], sep = sep.responses)
    resp.water <- cbind(WU.dat, as.data.frame(resp.water))
  } else
    resp.water <- WU.dat
  resp.water <- resp.water[do.call(order, resp.water),]
  return(resp.water)
}
