#Functions for calculating derived responses
#Function to get the values for specified times and, optionally, a column with the times factor 
"getTimesSubset" <- function(responses, times.factor = "Days", data, 
                       which.times, suffix = NULL, include.times.factor = FALSE,
                       include.individuals = FALSE, individuals = "Snapshot.ID.Tag")
{ 
  n <- dim(data)[1]
  subset <- rep(FALSE, n)
  for (day in which.times)
  { 
    this.times <- data[times.factor] == day
    subset <- subset | this.times
  }
  data.sub <- data[subset, ]
  data.sub <- data.sub[do.call(order, data.sub), ]
  if (!is.null(suffix))
    new.responses <- unlist(lapply(responses, 
                                   function(name, suffix){paste(name,suffix,sep=".")}, 
                                   suffix=suffix))
  else
    new.responses <- responses
  if (include.individuals)
  {
    responses <- c(responses,individuals)
    new.responses <- c(new.responses,individuals)
  }
  if (include.times.factor)
  {
    responses <- c(responses,times.factor)
    if (!is.null(suffix))
      new.responses <- c(new.responses, paste(times.factor,suffix,sep="."))
    else
      new.responses <- c(new.responses,times.factor)
  }
  data.sub <- data.sub[responses]
  names(data.sub) <- new.responses 
  return(data.sub)
}

#Function to form the growth rates over an interval for a set of responses
"intervalGRdiff" <- function(responses, individuals = "Snapshot.ID.Tag", 
                             which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                             times.factor = "Days", start.times, end.times, 
                             suffix.interval, data)
{ options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  
  #Check that individuals and times.factor are in data
  if (!all(c(individuals,times.factor) %in% names(data)))
    stop("Indivduals and/or times.factor are not in data")

  interval.resp <- cbind(getTimesSubset(responses, data = data, which.times = start.times, 
                                        times.factor = times.factor, 
                                        include.times.factor = TRUE, 
                                        suffix = "start"),
                         getTimesSubset(responses, data = data, which.times = end.times, 
                                        times.factor = times.factor, 
                                        include.times.factor = TRUE, 
                                        suffix = "end"))
  times.fac.start <- paste(times.factor,"start",sep=".")
  times.fac.end <- paste(times.factor,"end",sep=".")
  interval.resp[times.fac.start] <- as.numfac(interval.resp[[times.fac.start]])
  interval.resp[times.fac.end] <- as.numfac(interval.resp[[times.fac.end]])
  growth.rates <- lapply(responses, 
                         function(name, interval.resp, start.times, end.times, 
                                  which.rates, times.factor = "Days")
                         { #check which.rates
                           if (any(is.na(pmatch(which.rates, c("AGR","PGR","RGR")))))
                             stop("which.rates has at least one illegal value")
                           start.name <- paste(name,"start",sep=".")
                           start.times <- paste(times.factor,"start",sep=".")
                           end.name <- paste(name,"end",sep=".")
                           end.times <- paste(times.factor,"end",sep=".")
                           rates <- vector("list", length = length(which.rates))
                           k <- 0
                           if ("AGR" %in% which.rates)
                           { k <- k + 1
                             AGR <- (interval.resp[end.name] - interval.resp[start.name]) / 
                                          (interval.resp[end.times] - interval.resp[start.times])
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
                               (interval.resp[end.times] - interval.resp[start.times])
                             rates[k] <- RGR
                             if (is.null(suffices.rates))
                               names(rates)[k] <- paste(name,"RGR",suffix.interval,sep=".")
                             else
                               names(rates)[k] <- paste(name, suffices.rates[match("RGR",opt)],
                                                        suffix.interval, sep=".")
                           }
                           if ("PGR" %in% which.rates)
                           { k <- k +1
                             if ("RGR" %in% which.rates)
                               PGR <- exp(RGR)
                             else
                             { PGR <- (log(interval.resp[end.name]) - log(interval.resp[start.name])) / 
                                 (interval.resp[end.times] - interval.resp[start.times])
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
                         start.times = start.times, end.times = end.times, 
                         times.factor=times.factor, which.rates = opt)
  growth.rates <- as.data.frame(growth.rates)
  growth.rates <- cbind(growth.rates, 
                        getTimesSubset(individuals, data = data, times.factor = times.factor, 
                                       which.times = start.times, suffix = NULL))
  return(growth.rates)
}

#Function to form the growth rates over an interval for a set of responses by averaging growth rates
"intervalGRaverage" <- function(responses, individuals = "Snapshot.ID.Tag", 
                                which.rates = c("AGR","RGR"), suffices.rates=c("AGR","RGR"), 
                                start.time, end.time, times.factor = "Days", suffix.interval, 
                                data, sep=".", na.rm=TRUE)
{  options <- c("AGR","RGR")
   opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
   if (length(opt) != length(suffices.rates))
     stop("The length of of which.rates and suffices.rates should be equal")
   #Check that required growth rates are in data
   response.grates <- unlist(lapply(suffices.rates,
                                    function(grate, responses)
                                    { response.grates <- paste(responses,grate,sep=".")
                                      if (!all(response.grates %in% names(data)))
                                        stop("Growth rates for at least some responses are not in data")
                                      return(response.grates)
                                    },
                                    responses=responses))
   #Check that individuals and times.factor are in data
   if (!all(c(individuals,times.factor) %in% names(data)))
     stop("Indivduals and/or times.factor are not in data")
   #Get data for the times 
   times.vals <- unique(as.numfac(data[[times.factor]]))
   times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
   interval.resp <- getTimesSubset(response.grates, data = data, which.times = times.vals, 
                                   times.factor = times.factor, include.times.factor = TRUE)
   interval.resp <- cbind(getTimesSubset(individuals, data = data, which.times = times.vals, 
                                         times.factor = times.factor, 
                                         include.times.factor = FALSE),
                          interval.resp)
   interval.resp[times.factor] <- as.numfac(interval.resp[[times.factor]])
   #calculate the weights
   interval.resp <- split(interval.resp, as.list(interval.resp[individuals]), sep=sep)
   interval.resp <- lapply(interval.resp, 
                           function(data, times.factor = "Days")
                           { n <- nrow(data)
                             n1 <- n - 1
                             data[2:n, times.factor] <- data[2:n, times.factor] - 
                                                                   data[1:n1, times.factor]
                             data[1, times.factor] <- data[2, times.factor]
                             data[2:n1, times.factor] <- (data[1:(n1-1), times.factor] + 
                                                            data[2:n1, times.factor])/2
                             return(data)
                            }, times.factor = times.factor)
   interval.resp <- do.call(rbind, interval.resp)

   #Calculate growth rates by averaging within an interval
   avgrowth.rates <- lapply(responses, 
                            function(name, interval.resp, which.rates = c("AGR","RGR"))
                            { #check which.rates
                              if (any(is.na(pmatch(which.rates, c("AGR","RGR")))))
                                stop("which.rates has at least one illegal value")
                              k <- 0
                              if ("AGR" %in% which.rates)
                              { k <- k + 1
                                j <- match("AGR", which.rates)
                                name.AGR <- paste(name,suffices.rates[j],sep=".")
                                rates <- splitValueCalculate(name.AGR, weights=times.factor, 
                                                      individuals = individuals, 
                                                      FUN = "weighted.mean", data=interval.resp, 
                                                      na.rm=na.rm, sep=sep)
                                names(rates)[length(rates)] <- paste(name.AGR, suffix.interval,
                                                                     sep=".")
                              }
                              if ("RGR" %in% which.rates)
                              { k <- k + 1
                                j <- match("RGR", which.rates)
                                name.RGR <- paste(name,suffices.rates[j],sep=".")
                                interval.resp[name.RGR] <- log(interval.resp[name.RGR])
                                if (k==1)
                                  rates <- splitValueCalculate(name.RGR, weights=times.factor, 
                                                               individuals = individuals, 
                                                               FUN = "weighted.mean", 
                                                               data=interval.resp, 
                                                               na.rm=na.rm, sep=sep)
                                else
                                  rates <- merge(rates, 
                                                 splitValueCalculate(name.RGR, 
                                                                     weights=times.factor, 
                                                                     individuals = individuals, 
                                                                     FUN = "weighted.mean", 
                                                                     data=interval.resp, 
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
"intervalValueCalculate" <- function(response, weights=NULL, individuals = "Snapshot.ID.Tag", 
                                     FUN = "max", which.obs = FALSE, which.values = NULL, 
                                     start.time=NULL, end.time=NULL, times.factor = "Days", 
                                     suffix.interval=NULL, data, sep=".", na.rm=TRUE, ...)
{  
  #Trap which.levels and give message to replace it with which.values is not set
  tempcall <- list(...)
  if (length(tempcall) && "which.levels" %in% names(tempcall))
    stop("replace which.levels with which.values")
  
  #Check that response is in data
  if (!all(c(response,individuals,times.factor) %in% names(data)))
    stop("Some of the columns for response, indivduals and times.factor are not in data")
  #Get data for the times
  if (all(is.null(c(start.time, end.time))))
    interval.resp <- data[c(individuals,response,times.factor)]
  else
  { 
    times.vals <- unique(as.numfac(data[[times.factor]]))
    if (is.null(start.time))
      times.vals <- times.vals[times.vals <= end.time]
    else
      if (is.null(end.time))
        times.vals <- times.vals[times.vals >= start.time]
      else
        times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
      interval.resp <- getTimesSubset(response, data = data, which.times = times.vals, 
                                      times.factor = times.factor, include.times.factor = TRUE)
      interval.resp <- cbind(getTimesSubset(individuals, data = data, which.times = times.vals, 
                                            times.factor = times.factor, 
                                            include.times.factor = FALSE),
                             interval.resp)
  }
  
  #Calculate a value within an interval for each individual
  val.dat <- splitValueCalculate(response=response, weights=weights, individuals = individuals, 
                                 FUN = FUN, which.obs = which.obs, which.values = which.values, 
                                 data = interval.resp, na.rm=na.rm, sep=sep, ...)
  if (!is.null(suffix.interval))
  {
    names(val.dat)[match(paste(response,FUN,sep="."), names(val.dat))] <- 
      paste(response,FUN,suffix.interval,sep=".")
    if (which.obs)
      names(val.dat)[match(paste(response,FUN,"obs",sep="."), names(val.dat))] <- 
        paste(response,FUN,"obs",suffix.interval,sep=".")
    if (!is.null(which.values))
      names(val.dat)[match(paste(response,FUN,which.values,sep="."), names(val.dat))] <- 
        paste(response,FUN,which.values,suffix.interval,sep=".")
  }
  
  return(val.dat)
}

#Functions to calculate a single-valued function, including the observation has the value of the function

"which.funct.value" <- function(x, FUNCT = NULL, na.rm = TRUE, ...)
{ 
  funct <- get(FUNCT)
  funct <- match.fun(funct)
#  k <- which(x==funct(x, ...))
  adiff <- abs(x - funct(x, na.rm = na.rm, ...))
  k <- which(adiff == min(adiff, na.rm = na.rm))
  if (length(k) == 0)
    k <- NA
  else if (length(k) > 1)
    k <- k[1]
  return(k)
}
#Functions to calculate statistics for a profile

#Function to calculate the sample range
#"sample.range" <- function(x, na.rm = FALSE){z <- diff(range(x, na.rm = na.rm))}

#Function to calculate the root mean square
#"rms" <- function(x, na.rm = FALSE){z <- sqrt(sum(x*x, na.rm = na.rm)/length(x))}

#Function to calculate the sum over an interval between start.value and end.value
"interval.sum" <- function(x, start.value=1, end.value=1, na.rm = FALSE)
  # Get the sum of x for which subset.var has values between start.values and end.values
{ z <- sum(x[start.value:end.value], na.rm = na.rm)
}

#Function to calculate the sum over an interval between start.value and end.value
"interval.wtmean" <- function(x, w, start.value=1, end.value=1, na.rm = FALSE)
  # Get the sum of x for which subset.var has values between start.values and end.values
{ z <- sum(w[start.value:end.value] * x[start.value:end.value], na.rm = na.rm) /
           sum(w[start.value:end.value])
}

"splitValueCalculate" <- function(response, weights=NULL, individuals = "Snapshot.ID.Tag", 
                                  FUN = "max", which.obs = FALSE, which.values = NULL, 
                                  data, na.rm=TRUE, sep=".", ...)
  #a function to compute a FUN from the response for each individual
  #response is a character string giving the name of the response in data
  #individuals is a character vector giving the factors that index the individuals 
  #   for each of which a single value of funct is obtained from their observations
  #... allows for optional arguments to FUN
{ 
  #Trap which.levels and give message to replace it with which.values is not set
  tempcall <- list(...)
  if (length(tempcall) && "which.levels" %in% names(tempcall))
    stop("replace which.levels with which.values")
 
  funct <- get(FUN)
  funct <- match.fun(funct)
  #Check that response and individuals are in data
  if (!all(c(response, individuals) %in% names(data)))
    stop("response and/or indivduals are not in data")
  
  nfac <- length(individuals)
  if (!nfac) 
    stop("'individuals' is of length zero")
  kresp <- match(response, names(data))
  #Form data frame of values
  data <- split(data, as.list(data[individuals]))
  if (is.null(weights))
    val.dat <- lapply(data, 
                       function(data, response, FUNC, na.rm, ...)
                       { 
                         vals <- FUNC(x=data[[response]], na.rm=na.rm, ...) 
                         return(vals)
                       },
                       response=response, FUNC=funct, na.rm=na.rm, ...)
  else
    val.dat <- lapply(data, 
                       function(data, response, weights, FUNC, na.rm, ...)
                       { 
                         vals <- FUNC(x=data[[response]], w= data[[weights]], na.rm=na.rm, ...) 
                         return(vals)
                       },
                       response=response, weights=weights, FUNC=funct, na.rm=na.rm, ...)
  val.dat <- as.data.frame(do.call(rbind, val.dat))
  names(val.dat) <- paste(response,FUN,sep=".")
  val.dat[[1]][is.infinite(val.dat[[1]])] <- NA
  indices <- rownames(val.dat)
  indices <- strsplit(indices, split=sep, fixed=TRUE)
  for (fac in 1:length(individuals))
  { 
    val.dat[[individuals[fac]]] <- unlist(lapply(indices, 
                                              function(x, fac)
                                                { x[fac]}, 
                                              fac))
    if (is.factor(data[[individuals[fac]]]))
      val.dat[[individuals[fac]]] <- factor(val.dat[[individuals[fac]]])
    else
      if (is.numeric(data[[individuals[fac]]]))
        val.dat[[individuals[fac]]] <- as.numeric(val.dat[[individuals[fac]]])
  }
  val.dat <- val.dat[, c(2:length(val.dat),1)]
  
  #Get which observation is equal to each returned function value, if required
  if (which.obs | !is.null(which.values))
  { 
    kresp.val <- length(val.dat)
    resp.val <- names(val.dat)[kresp.val]
    which.dat <- lapply(data, 
                        function(x, response, FUNCT = NULL, na.rm = TRUE, 
                                 which.obs, which.values, ...)
                        { 
                          #Find which observation number corresponds to the value of FUNCT
                          w <- which.funct.value(x[[response]], FUNCT = FUNCT, na.rm = na.rm, ...)
                          #Match observation numbers with the corresponding values of the factor/numeric which.values
                          if (!is.null(which.values))
                          {
                            val <- x[[which.values]][w]
                            if (which.obs)
                              w <- list(w,val)
                            else
                              w <- list(val)
                          } else
                            w <- list(w)
                          return(w)
                        },
                        response=response, FUNCT = FUN, na.rm = na.rm, 
                        which.obs = which.obs, which.values = which.values, ...)
    if (which.obs && !is.null(which.values))
    {
      which.dat <- do.call(rbind, 
                           lapply(which.dat, 
                                  function(x)
                                  {
                                    x <- data.frame(x)
                                    names(x) <- c("V1","V2")
                                    return(x)
                                  }))
      names(which.dat) <- paste(resp.val, c("obs", which.values), sep=".")
      ntab <- length(which.dat)
      which.dat[[ntab-1]][is.infinite(which.dat[[ntab-1]])] <- NA
      which.dat[[ntab]][is.infinite(which.dat[[ntab]])] <- NA
    } else
    {
      which.dat <- as.data.frame(unlist(which.dat))
      if (!is.null(which.values))
        resp.which <- paste(resp.val,which.values,sep=".") 
      else
        resp.which <- paste(resp.val,"obs",sep=".")
      names(which.dat) <-  resp.which
      ntab <- length(which.dat)
      which.dat[[ntab]][is.infinite(which.dat[[ntab]])] <- NA
      kresp.val <- kresp.val + 1
    }
    rownames(val.dat) <- rownames(which.dat) <- NULL
    val.dat <- cbind(val.dat,which.dat)
  }
  
  #Put data frame into standard order
  val.dat <- val.dat[do.call(order, val.dat), ]
  
  return(val.dat)
}

# Function to calculate water use indices (WUI) over an interval for a set of responses
"intervalWUI" <- function(responses, water.use = "Water.Use", individuals = "Snapshot.ID.Tag", 
                          times.factor = "Days", start.times, end.times, suffix.interval = NULL, 
                          data, include.total.water = FALSE, na.rm = FALSE)
{   #Check that response and individuals are in data
  if (!all(c(water.use, individuals, times.factor) %in% names(data)))
    stop("One or more of water use, response and indivduals are not in data")
  
  #get the water use
  sum.dat <- splitValueCalculate(water.use, individuals = individuals, 
                                 FUN = "sum", na.rm = na.rm, 
                                 data = subset(data, 
                                               as.numfac(eval(parse(text=times.factor))) >= 
                                                 max(start.times)+1 & 
                                               as.numfac(eval(parse(text=times.factor))) <= 
                                                 max(end.times)))
  if (!is.null(suffix.interval))
  { 
    water.name <- paste(water.use,"Total",suffix.interval,sep=".")
    wui.name <- paste("WUI",suffix.interval,sep=".")
  }
  else
  { 
    water.name <- paste(water.use,"Total",sep=".")
    wui.name <- paste("WUI",sep=".")
  }
  names(sum.dat)[match(paste(water.use,"sum",sep="."), names(sum.dat))] <- water.name
  
  #Get the values of the responses for the start and end of the time interval
  interval.resp <- cbind(getTimesSubset(individuals, data = data, times.factor = times.factor, 
                                        which.times = start.times, suffix = NULL),
                         getTimesSubset(responses, data = data, times.factor = times.factor, 
                                        which.times = start.times, include.times.factor = TRUE, 
                                        suffix = "start"),
                         getTimesSubset(responses, data = data, times.factor = times.factor, 
                                        which.times = end.times, include.times.factor = TRUE, 
                                        suffix = "end"))
  interval.resp <- merge(interval.resp, sum.dat, by = individuals)
  interval.resp <- interval.resp[ do.call(order, interval.resp), ]
  water.index <- lapply(responses, 
                        function(name, interval.resp, start.times, end.times, water.name, 
                                 include.total.water = FALSE)
                        { 
                          start.name <- paste(name,"start",sep=".")
                          end.name <- paste(name,"end",sep=".")
                          if (include.total.water)
                          { 
                            rates <- vector("list", length = 2)
                            rates[[1]] <- interval.resp[end.name][,1] - 
                                                        interval.resp[start.name][,1]
                            rates[[2]] <- WUI(rates[[1]], interval.resp[water.name][,1])
                            if (!is.null(suffix.interval))
                            { 
                              names(rates)[1] <- paste(name,"Total",suffix.interval,sep=".")
                              names(rates)[2] <- paste(name,"WUI",suffix.interval,sep=".")
                            } else
                            { 
                              names(rates)[1] <- paste(name,"Total",sep=".")
                              names(rates)[2] <- paste(name,"WUI",sep=".")
                            }
                          } else
                          { 
                            rates <- vector("list", length = 1)
                            rates[[1]] <- WUI((interval.resp[end.name][,1] - 
                                                 interval.resp[start.name][,1]),  
                                              interval.resp[water.name][,1])
                            if (!is.null(suffix.interval))
                              names(rates)[1] <- paste(name,"WUI",suffix.interval,sep=".")
                            else
                              names(rates)[1] <- paste(name,"WUI",sep=".")
                          }
                          rates <- as.data.frame(rates)
                          return(rates)
                        },
                        interval.resp = interval.resp, 
                        start.times = start.times, end.times = end.times, 
                        water.name = water.name, include.total.water = include.total.water)
  water.index <- as.data.frame(water.index)
  if (include.total.water)
    water.index <- cbind(interval.resp[c(individuals, water.name)], water.index)
  else 
    water.index <- cbind(interval.resp[individuals], water.index)
  return(water.index)
}
