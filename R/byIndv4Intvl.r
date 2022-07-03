#Functions to produce a value in an interval for each interval
#Function to form the growth rates over an interval for a set of responses
"byIndv4Intvl_GRsDiff" <- function(data, responses, 
                                   individuals = "Snapshot.ID.Tag", times = "DAP", 
                                   which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                                   start.time, end.time, suffix.interval)
{ 
  options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  
  if (!all(sapply(list(start.time, end.time, suffix.interval), 
                  function(x) length(x) == 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length one")
  
  #Check that responses, individuals and times are in data
  checkNamesInData(c(responses, individuals, times), data = data)
  
  interval.resp <- cbind(getTimesSubset(data = data, responses = responses, 
                                        times = times, which.times = start.time, 
                                        include.times = TRUE, 
                                        suffix = "start"),
                         getTimesSubset(data = data, responses = responses, 
                                        times = times, which.times = end.time, 
                                        include.times = TRUE, 
                                        suffix = "end"))
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
                               names(rates)[k] <- paste(name,"AGR",suffix.interval,sep=".")
                             else
                               names(rates)[k] <- paste(name, suffices.rates[match("AGR",opt)],
                                                        suffix.interval, sep=".")
                           }
                           if ("RGR" %in% which.rates)
                           { k <- k + 1
                           RGR <- (log(interval.resp[end.name]) - log(interval.resp[start.name])) / 
                             (interval.resp[end.time] - interval.resp[start.time])
                           rates[k] <- RGR
                           if (is.null(suffices.rates))
                             names(rates)[k] <- paste(name,"RGR",suffix.interval,sep=".")
                           else
                             names(rates)[k] <- paste(name, suffices.rates[match("RGR",opt)],
                                                      suffix.interval, sep=".")
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
                               names(rates)[k] <- paste(name,"PGR",suffix.interval,sep=".")
                             else
                               names(rates)[k] <- paste(name, suffices.rates[match("PGR",opt)],
                                                        suffix.interval, sep=".")
                           }
                           rates <- as.data.frame(rates)
                           return(rates)
                         },
                         interval.resp = interval.resp, 
                         start.time = start.time, end.time = end.time, 
                         times=times, which.rates = opt)
  growth.rates <- as.data.frame(growth.rates)
  #use getTimesSubset to add the individuals to the growth.rates data.frame  
  growth.rates <- cbind(growth.rates, 
                        getTimesSubset(data = data, responses = individuals, times = times, 
                                       which.times = start.time, suffix = NULL))
  return(growth.rates)
}

#Function to form the growth rates over an interval for a set of responses by averaging growth rates
"byIndv4Intvl_GRsAvg" <- function(data, responses, 
                                  individuals = "Snapshot.ID.Tag", times = "DAP", 
                                  which.rates = c("AGR","RGR"), suffices.rates=c("AGR","RGR"), 
                                  start.time, end.time, suffix.interval, 
                                  sep=".", na.rm=TRUE)
{  
  #Check that responses, individuals and times are in data
  checkNamesInData(c(responses, individuals, times), data = data)

  if (!all(sapply(list(start.time, end.time, suffix.interval), 
                  function(x) length(x) == 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length one")
  
  options <- c("AGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  #Check that required growth rates are in data
  response.grates <- unlist(lapply(suffices.rates,
                                   function(grate, responses)
                                   { 
                                     response.grates <- paste(responses,grate,sep=".")
                                     if (!all(response.grates %in% names(data)))
                                       stop("The following growth rates are not in data: ",
                                            paste0(response.grates[!(response.grates %in% names(data))]))
                                     return(response.grates)
                                   },
                                   responses=responses))
  #Get data for the times 
  times.vals <- unique(convertTimes2numeric(data[[times]]))
  times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
  interval.resp <- getTimesSubset(data = data, responses = response.grates, 
                                  times = times, which.times = times.vals, 
                                  include.times = TRUE)
  interval.resp <- cbind(getTimesSubset(data = data, responses = individuals, 
                                        times = times, which.times = times.vals, 
                                        include.times = FALSE),
                         interval.resp)
  interval.resp[times] <- convertTimes2numeric(interval.resp[[times]])
  #calculate the weights
  interval.resp <- split(interval.resp, as.list(interval.resp[individuals]), sep=sep)
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
                             if (any(is.na(pmatch(which.rates, c("AGR","RGR")))))
                               stop("which.rates has at least one illegal value")
                             k <- 0
                             if ("AGR" %in% which.rates)
                             { 
                               k <- k + 1
                               j <- match("AGR", which.rates)
                               name.AGR <- paste(name,suffices.rates[j],sep=".")
                               rates <- byIndv_ValueCalc(data=interval.resp, response = name.AGR, 
                                                         individuals = individuals, 
                                                         FUN = "weighted.mean", weights=times, 
                                                         na.rm=na.rm, sep=sep)
                               names(rates)[length(rates)] <- paste(name.AGR, suffix.interval,
                                                                    sep=".")
                             }
                             if ("RGR" %in% which.rates)
                             { 
                               k <- k + 1
                               j <- match("RGR", which.rates)
                               name.RGR <- paste(name,suffices.rates[j],sep=".")
                               interval.resp[name.RGR] <- log(interval.resp[name.RGR])
                               if (k==1)
                                 rates <- byIndv_ValueCalc(data=interval.resp, response = name.RGR, 
                                                           individuals = individuals, 
                                                           FUN = "weighted.mean", 
                                                           weights=times, 
                                                           na.rm=na.rm, sep=sep)
                               else
                                 rates <- merge(rates, 
                                                byIndv_ValueCalc(data=interval.resp, response = name.RGR, 
                                                                 individuals = individuals, 
                                                                 FUN = "weighted.mean", 
                                                                 weights=times, 
                                                                 na.rm=na.rm, sep=sep),
                                 )
                               rates[length(rates)] <- exp(rates[length(rates)])
                               names(rates)[length(rates)] <- paste(name.RGR, suffix.interval,
                                                                    sep=".")
                             }
                             return(rates)
                           },
                           interval.resp = interval.resp, 
                           which.rates = opt)
  avgrowth.rates <- as.data.frame(avgrowth.rates)
  return(avgrowth.rates)
}

#Function to calculate a value for observations within an interval for a set of responses
"byIndv4Intvl_ValueCalc" <- function(data, response, 
                                     individuals = "Snapshot.ID.Tag", times = "DAP", 
                                     FUN = "max", which.obs = FALSE, which.values = NULL, 
                                     addFUN2name = TRUE, 
                                     start.time=NULL, end.time=NULL, suffix.interval=NULL, 
                                     sep=".", weights=NULL, na.rm=TRUE, ...)
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
                              addFUN2name = addFUN2name, weights=weights, na.rm=na.rm, sep=sep, ...)
  if (!is.null(suffix.interval))
  {
    if (addFUN2name)
    {   
      names(val.dat)[match(paste(response, FUN, sep="."), names(val.dat))] <- 
        paste(response, FUN, suffix.interval, sep=".")
      if (which.obs)
        names(val.dat)[match(paste(response,FUN, "obs", sep="."), names(val.dat))] <- 
          paste(response, FUN, "obs", suffix.interval, sep=".")
      if (!is.null(which.values))
        names(val.dat)[match(paste(response, FUN, which.values, sep="."), names(val.dat))] <- 
          paste(response, FUN, which.values, suffix.interval, sep=".")
    } else
    {
      names(val.dat)[match(response, names(val.dat))] <- 
        paste(response, suffix.interval, sep=".")
      if (which.obs)
        names(val.dat)[match(paste(response, "obs", sep="."), names(val.dat))] <- 
          paste(response, "obs", suffix.interval, sep=".")
      if (!is.null(which.values))
        names(val.dat)[match(paste(response, which.values, sep="."), names(val.dat))] <- 
          paste(response, which.values, suffix.interval, sep=".")
    }
  }
  
  return(val.dat)
}

# Function to calculate water use indices (WUI) over an interval for a set of responses
"byIndv4Intvl_WaterUse" <- function(data, water.use = "Water.Use", responses = NULL, 
                                    individuals = "Snapshot.ID.Tag", times = "DAP", 
                                    trait.types = c("WU", "WUR", "WUI"), 
                                    suffix.rate = "R", suffix.index = "I", 
                                    start.time, end.time, suffix.interval = NULL, 
                                    na.rm = FALSE)
{ 
  options <- c("WU", "WUR", "AGR", "WUI", "all")
  traits <- options[unlist(lapply(trait.types, check.arg.values, options=options))]
  if ("all" %in% traits)
    traits <- c("WU", "WUR", "AGR", "WUI")
  traits <- c("WU", "WUR", "AGR", "WUI")[c("WU", "WUR", "AGR", "WUI") %in% traits]
  names(traits) <- traits
  if (length(traits) == 1 && traits == "AGR")
    warning("no water use traits have specified - only AGR")

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
  trait.names <- c("WU", "WUR", "AGR", "WUI")
  names(trait.names) <- trait.names
  trait.names["WU"] <- water.use
  trait.names["WUR"] <- paste0(water.use, suffix.rate)
  trait.names["WUI"] <- paste0(water.use, suffix.index)
  trait.names <- trait.names[traits]
  trait.names <- paste(trait.names, suffix.interval, 
                       sep = ifelse(is.null(suffix.interval), "", "."))
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
  
  #Get the AGR and WUI traits, as requested
  if (any(c("AGR", "WUI") %in% traits))
  {  
    resp.water <- lapply(responses, 
                         function(name, interval.resp, start.time, end.time, trait.names, 
                                  traits)
                         { 
                           start.name <- paste(name,"start",sep=".")
                           end.name <- paste(name,"end",sep=".")
                           #Calculate the water traits
                           rates <- list()
                           if (any(c("AGR", "WUI") %in% traits))
                           {
                             rates <- c(rates, list((interval.resp[end.name][,1] - 
                                                       interval.resp[start.name][,1])))
                             if ("AGR" %in% traits)
                               names(rates)[length(rates)] <- paste(name, trait.names["AGR"], sep = ".")
                             if ("WUI" %in% traits)
                             {
                               if ("AGR" %in% traits)
                                 rates <- c(rates, rates[length(rates)])
                               kWUI <- length(rates)
                               rates[[kWUI]] <- WUI(rates[[kWUI]], 
                                                    interval.resp["WU"][,1])
                               names(rates)[kWUI] <- paste(name, trait.names["WUI"], sep = ".")
                             }
                           }
                           rates <- as.data.frame(rates)
                           return(rates)
                         },
                         interval.resp = interval.resp, 
                         start.time = start.time, end.time = end.time, 
                         trait.names = trait.names, traits = traits)
    resp.water <- cbind(WU.dat, as.data.frame(resp.water))
    } else
      resp.water <- WU.dat
  return(resp.water)
}
