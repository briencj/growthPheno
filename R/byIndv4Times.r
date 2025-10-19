#Functions to produce values for all times for each individual
"byIndv4Times_GRsDiff" <- function(data, responses, 
                                   individuals = "Snapshot.ID.Tag", times = "DAP", 
                                   which.rates = c("AGR","PGR","RGR"), 
                                   suffices.rates=NULL, sep.rates = ".", 
                                   avail.times.diffs = FALSE, ntimes2span = 2)
{ 
  options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  
  #Check that responses, individuals and times are in data
  vars <- c(individuals, times, responses)
  checkNamesInData(vars, data = data)

  #Deal with time.diffs
  times.diffs <- paste(times, "diffs", sep=".")
  if (avail.times.diffs)
  { 
    if (times.diffs %in% names(data))
      vars <- c(vars, times.diffs)
    else
      stop("The column ", times.diffs, ", expected to contain the times differences, is not in data")
  } else #Check that all DAPs available for all plants
  {
    DAP.reps <- table(data[times])
   if (all(length(unique(data[individuals][[1]])) != DAP.reps[1]))
     warning("Not all individuals have an entry for the first of the times, as is assumed in calculating \ntime differences; ",
             "either add the absent times with NA for the observed response values or \n",
             "add your own column of differences to the data and use the avail.time.diffs argument.")
  }
  
  tmp <- data[vars]
  tmp <- tmp[do.call(order, tmp), ]
  lag <- ntimes2span - 1
  xTime <- convertTimes2numeric(tmp[[times]])
  miss.days <- sort(unique(xTime))[1:lag]
  
  #time.diffs are available - check them
  if (avail.times.diffs && 
      !all(is.na(data[[times.diffs]][xTime %in% miss.days])))
    stop("The times.diffs column in data does not have the appropriate missing values for the ", 
         "inital times of each individual\n",
         "Set avail.times.diffs to FALSE to have 'byIndv4Intvl_GRsDiff' calculate them")
  
  if (any(is.na(data[[times]])))
    warning(paste("Some values of ",times,
                  " are missing, which can result in merge producing a large data.frame", 
                  sep = ""))
  if (any(unlist(lapply(as.list(data[individuals]), 
                        function(f)
                          any(is.na(f))))))
    warning(paste("Some values of the factors in individuals are missing, ",
                  "which can result in merge producing a large data.frame", sep = ""))
  
  #Form time differences in a way that every first time is the same Time
  # - setting first time point to missing results in the growth rates also being NA
  if (!avail.times.diffs)
  { 
    #nspan (t)  2   3    4  5    6  7
    #nmiss      1   2    3  4    5  6 = lag
    #posn1      2   2    3  3    4  4 = ceiling((t+1)/2)
    #(t+1)/2  1.5   2  2.5  3  3.5  4
    #NA move    0   1    1  2    2  3 = nmiss - (posn1 -1) = floor((t+1) /2) - 1
    
    tmp[times.diffs] <- calcLagged(xTime, operation ="-", lag = lag)
    #Make sure that first times.diffs are NA
    if (!all(is.na(tmp[[times.diffs]][xTime %in% miss.days])))
      tmp[[times.diffs]][xTime %in% miss.days] <- NA
    rownames(tmp) <- NULL
  }
  
  responses.GRs <- c()
  #Form AGR (because time.diffs is NA for first time, so will the growth rates)
  if ("AGR" %in% opt)
  { 
    if (is.null(suffices.rates))
      responses.GR <- paste(responses, "AGR", sep=sep.rates)
    else
      responses.GR <- paste(responses, suffices.rates[match("AGR",opt)], sep=sep.rates)
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                              FUN = AGRdiff, 
                                              time.diffs = tmp[[times.diffs]], lag = lag))
    responses.GRs <- c(responses.GRs, responses.GR)
  }
  
  #Form PGR (because time.diffs is NA for first time, so will the growth rates)
  if ("PGR" %in% opt)
  { 
    if (is.null(suffices.rates))
      responses.GR <- paste(responses, "PGR", sep=sep.rates)
    else
      responses.GR <- paste(responses, suffices.rates[match("PGR",opt)], sep=sep.rates)
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                              FUN = PGR, 
                                              time.diffs = tmp[[times.diffs]], lag = lag))
    responses.GRs <- c(responses.GRs, responses.GR)
  }
  
  #Form RGR (because time.diffs is NA for first time, so will the growth rates)
  if ("RGR" %in% opt)
  { 
    if (is.null(suffices.rates))
      responses.GR <- paste(responses, "RGR", sep=sep.rates)
    else
      responses.GR <- paste(responses, suffices.rates[match("RGR",opt)], sep=sep.rates)
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                              FUN = RGRdiff, 
                                              time.diffs = tmp[[times.diffs]], lag = lag))
    responses.GRs <- c(responses.GRs, responses.GR)
  }
  
  #Reposition the GRs, if necessary
  nmove <- floor((ntimes2span+1) /2) - 1
  if (nmove > 0)
  {
    if (!avail.times.diffs)
      cols2move <- c(times.diffs, responses.GRs)
    else
      cols2move <- responses.GRs
    tmp <- split(tmp, f = as.list(tmp[individuals]))
    tmp <- lapply(tmp,
                  function(x, cols2move, nmove)
                  {
                    x[cols2move] <- x[c((nmove+1):nrow(x), 1:nmove), cols2move]
                    return(x)
                  },
                  cols2move = cols2move, nmove = nmove)
    tmp <- do.call(rbind, tmp)
  }
  
  #Remove NAs in individuals and time.factor in tmp
  if (any(is.na(tmp[[times]])))
    tmp <- tmp[!is.na(tmp[[times]]), ]
  if (any(unlist(lapply(as.list(tmp[individuals]), 
                        function(f)
                          any(is.na(f))))))
    for (f in individuals)
      tmp <- tmp[!is.na(tmp[f]),]
  tmp <- tmp[,-match(responses,names(tmp))]
  #Keep times.diffs in data if they were used
  if (avail.times.diffs)
    tmp <- tmp[,-match(times.diffs, names(tmp))]
  #Remove unused times.diffs from data so time.diffs used are added from tmp
  if (!avail.times.diffs && times.diffs %in% names(data))
    data <- data[,-match(times.diffs, names(data))]
  if (any(responses.GR %in% names(data)))
    data <- data[,-match(responses.GRs, names(data))]
  data <- left_join(data, tmp, by = c(individuals, times))
  data  <- data[do.call(order, data),]
  return(data)
}



#Fit splines to smooth the longitudinal trends in a set of individuals for a response
#Specify responses to be smoothed and then loop over the individuals
"indvSplines" <- function(data, response, response.smoothed, 
                          individuals, times, 
                          smethod, stype, df, lambda, npspline.segments, 
                          correctBoundaries, 
                          rates, suffices.rates, 
                          extra.derivs, suffices.extra.derivs, 
                          #deriv, suffices.deriv, extra.rate, 
                          na.x.action, na.y.action, sep.levels, ...)
{
  #Split data frame by each combination of the individuals factors
  old.names <- names(data)
  tmp <- split(data, as.list(data[individuals]), sep=sep.levels)
  #Fit splines for each combination of the individuals factors
  tmp <- lapply(tmp, smoothSpline, 
                response = response, response.smoothed = response.smoothed, 
                x = times, smoothing.method = smethod, spline.type = stype, 
                df=df, lambda = lambda, npspline.segments = npspline.segments, 
                correctBoundaries = correctBoundaries, 
                rates = rates, suffices.rates = suffices.rates,
                extra.derivs = extra.derivs, 
                suffices.extra.derivs = suffices.extra.derivs, 
#                suffices.deriv = suffices.deriv, 
#                extra.rate = extra.rate, 
                na.x.action = na.x.action, na.y.action = na.y.action, ...)
  tmp <- lapply(tmp, function(dat) dat$predictions) #extract predictions
  tmp <- do.call(rbind, tmp)
  ncols <- ncol(tmp)
  indices <- rownames(tmp)
  indices <- strsplit(indices, split=sep.levels, fixed=TRUE)
  for (fac in 1:length(individuals))
  { 
    tmp[[individuals[fac]]] <- unlist(lapply(indices, 
                                             function(x, fac)
                                             { x[fac]}, 
                                             fac))
    if (is.factor(data[[individuals[fac]]]))
      tmp[[individuals[fac]]] <- factor(tmp[[individuals[fac]]])
    else
      if (is.numeric(data[[individuals[fac]]]))
        tmp[[individuals[fac]]] <- as.numeric(tmp[[individuals[fac]]])
  }
  tmp <- tmp[, c((ncols+1):length(tmp),1:ncols)]
  #Remove any pre-existing smoothed cols in data
  tmp.smooth <- names(tmp)[-match(c(individuals,times), names(tmp))]
  tmp.smooth <- na.omit(match(tmp.smooth, names(data)))
  if (length(tmp.smooth) > 0)
    data <- data[ ,-tmp.smooth]
  tmp <- tmp[!is.na(tmp[[times]]), ]
  data <- merge(data, tmp, all.x = TRUE, sort=FALSE)
  #Rearrange columns so original column are followed by new columns
  new.names <- names(data)
  new.names <- new.names[-match(old.names, new.names)]
  data <- data[c(old.names, new.names)]
  rownames(data) <- NULL
  
  return(data)
}

#The function that controls the fitting
"byIndv4Times_SplinesGRs" <- function(data, response, response.smoothed = NULL, 
                                      individuals = "Snapshot.ID.Tag", times,  
                                      smoothing.method = "direct", smoothing.segments = NULL, 
                                      spline.type = "NCSS", df=NULL, lambda = NULL, 
                                      npspline.segments = NULL, 
                                      correctBoundaries = FALSE,
                                      rates.method = "differences", 
                                      which.rates = c("AGR","RGR"), 
                                      suffices.rates = NULL, sep.rates = ".",
                                      avail.times.diffs = FALSE, ntimes2span = 2, 
                                      extra.derivs = NULL, suffices.extra.derivs=NULL, 
                                      sep.levels=".", 
                                      na.x.action="exclude", na.y.action = "trimx", ...)
{ 

  impArgs <- match.call()
  if ("na.rm" %in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  if ("smoothing.scale" %in% names(impArgs))
    stop("smoothing.scale has been deprecated; use smoothing.method")
  
  options <- c("none", "differences","derivatives")
  ratemeth.opt <- options[check.arg.values(rates.method, options=options)]
  
  options <- c("AGR", "PGR", "RGR")
  grates <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (correctBoundaries && (ratemeth.opt != "none"))
    stop("Unable to correctBoundaries when rates.method is not none")
  if ((ratemeth.opt == "derivatives") && ("PGR" %in% grates))
    stop("PGR is not available when rates are based on derivatives")
  if (!is.null(suffices.rates) & length(grates) != length(suffices.rates))
    stop("The length of of rates and suffices.rates should be equal")
  if ((ratemeth.opt == "derivatives") && !is.null(extra.derivs) && (1 %in% extra.derivs))
    stop("when rates.method is derivatives, 1 should not be included in extra.derivs")
  if (is.null(suffices.rates))
    suffices.rates <- grates
  names(suffices.rates) <- grates
  
  #Check that response, individuals and times are in data
  vars <- c(individuals, times, response)
  checkNamesInData(vars, data = data)
  times.diffs <- paste(times, "diffs", sep=".")
  if (avail.times.diffs)
  { 
    if (times.diffs %in% names(data))
      vars <- c(vars, times.diffs)
    else
      stop("The column ", times.diffs, ", expected to contain times differences, is not in data")
  }
  
  smethods <- c("direct", "logarithmic")
  smethod <- smethods[check.arg.values(smoothing.method, options=smethods)]
  stype <- c("NCSS", "PS")
  stype <- stype[check.arg.values(spline.type, options=stype)]
  if (is.null(lambda) && stype == "PS")
    stop("Must specify lambda for spline.type set to PS")
  if (!is.null(df) && !is.null(lambda))
    stop("One of df and lambda must be NULL")
  
  if (is.null(response.smoothed))
    response.smoothed <- paste0("s", response)
  
  #Check npspline.segments
  if (length(npspline.segments) > 1)
  { 
    if (is.null(smoothing.segments))
      stop("npspline.segments must be of length one when smoothing.segments is NULL")
    else
    {
      if (length(npspline.segments) != length(smoothing.segments))
        stop("the number of values of npspline.segments should be one or ",
             "equal to the number of segments in a segmented spline fit")
    }
    if(!all(diff(unlist(smoothing.segments)) > 0))
      stop("the smoothing.segments are not a set of non-overlapping, successive intervals")
  }
  
  #Check that for overlapping smoothing segments and give an error or a warning depending on get.rates
  if (!all(diff(unlist(smoothing.segments)) > 0))
  { 
    if (!("none" %in% ratemeth.opt))
      stop("rates.method must be `none` when times values occur in more than one smoothing segment")
    else
      warning("The values for some times occur in multiple smoothing segments and so some individuals ",
              "will have multiple rows in the returned data.frame, one for each segment in which the ",
              "times occur.")
  }
  
  #Determine rates and suffices rates for spline fitting
  splrates <- NULL
  splsuffices <- NULL
  if (ratemeth.opt == "derivatives")
  {
    splrates <- grates
    splsuffices <- suffices.rates
  }

  #Fit the splines to the individuals
  tmp <- data
  tmp[times] <- convertTimes2numeric(tmp[[times]])
  if (is.allnull(smoothing.segments))
  { 
    smth <- indvSplines(data = tmp,
                        response=response, response.smoothed=response.smoothed, 
                        individuals = individuals, times = times, 
                        smethod = smethod, stype = stype, df=df, 
                        lambda = lambda, npspline.segments = npspline.segments, 
                        correctBoundaries = correctBoundaries, 
                        rates = splrates, 
                        suffices.rates = splsuffices, sep.rates = sep.rates,
                        extra.derivs = extra.derivs, 
                        suffices.extra.derivs = suffices.extra.derivs, 
#                        deriv = derivs, suffices.deriv = suffices.derivs, 
#                        extra.rate = extra.rate, 
                        na.x.action = na.x.action, na.y.action = na.y.action, 
                        sep.levels = sep.levels, ...)
    if (rates.method == "differences")
    {
      smth <- byIndv4Times_GRsDiff(data = smth, response.smoothed, 
                                   individuals=individuals,times=times, 
                                   which.rates = grates, sep.rates = sep.rates, 
                                   ntimes2span = ntimes2span, 
                                   avail.times.diffs = avail.times.diffs)
    }
    smth[times] <- convertTimesExnumeric(smth[[times]], data[[times]])
  } else    
  {
    knseg <- npspline.segments[1]
    smth <- data.frame()
    for (k in 1:length(smoothing.segments))
    {
      segm <- smoothing.segments[[k]]
      subdat <- tmp[(tmp[times] >= segm[1]) & (tmp[times] <= segm[2]),] 
      if (length(npspline.segments) > 1) knseg <- npspline.segments[k]
      subdat<- indvSplines(data = subdat, 
                           response=response, response.smoothed=response.smoothed, 
                           individuals = individuals, times = times, 
                           smethod = smoothing.method, stype = spline.type, 
                           lambda = lambda, df=df, npspline.segments = npspline.segments, 
                           correctBoundaries = correctBoundaries, 
                           rates = splrates, suffices.rates = splsuffices,
                           extra.derivs = extra.derivs, 
                           suffices.extra.derivs = suffices.extra.derivs, 
#                           deriv = derivs, suffices.deriv = suffices.derivs, 
#                           extra.rate = extra.rate, 
                           na.x.action=na.x.action, na.y.action=na.y.action, 
                           sep.levels = sep.levels, ...)
      if (rates.method == "differences" && ntimes2span != 2)
        smth <- byIndv4Times_GRsDiff(data = smth, response.smoothed, 
                                     individuals=individuals,times=times, 
                                     which.rates = grates, sep.rates = sep.rates, 
                                     ntimes2span = ntimes2span, 
                                     avail.times.diffs = avail.times.diffs)
      subdat[times] <- convertTimesExnumeric(subdat[[times]], data[[times]])
      smth <- rbind(smth, subdat)
    }
    smth <- smth[do.call(order, smth), ]
    if (rates.method == "differences" && ntimes2span == 2)
      smth <- byIndv4Times_GRsDiff(data = smth, response.smoothed, 
                                   individuals = individuals, times = times, 
                                   which.rates = grates, sep.rates = sep.rates, 
                                   ntimes2span = ntimes2span, 
                                   avail.times.diffs = avail.times.diffs)
  }
  
  return(smth)
}


"byIndv4Times_WaterUse" <- function(data, weight.after = "Weight.After",
                                  weight.before = NULL, water.added = NULL, 
                                  responseAGR = NULL, 
                                  individuals = "Snapshot.ID.Tag", times = "DAP", 
                                  which.trait.types = c("WU","WUR","WUI"), 
                                  water.trait.names = c("WU", "WUR", "WUI"),  
                                  avail.times.diffs = FALSE, 
                                  ntimes2span = 2)
{ 
  options <-c("WU","WUR","WUI")
  opt <- options[unlist(lapply(which.trait.types, check.arg.values, options=options))]
  
  if (all(is.null(c(weight.before, water.added))) && length(c(weight.before, water.added)) != 1)
    stop("One and only one of weight.before and water.added must be NULL")
  
  if (length(water.trait.names) < length(which.trait.types))
    stop("The length of water.trait.names must be at least that of which.trait.types.")
  else
    water.trait.names <- water.trait.names[1:(length(which.trait.types))]

    #Check that weight.after, weight.before, water.added, individuals and times are in data
  vars <- c(individuals, times, weight.after, weight.before, water.added, responseAGR)
  checkNamesInData(vars, data = data)
  
  #Check for time.diffs
  times.diffs <- paste(times, "diffs", sep=".")
  if (avail.times.diffs)
  { 
    if (times.diffs %in% names(data))
      vars <- c(vars, times.diffs)
    else
      stop("The column ", times.diffs, ", expected to contain the times differences, is not in data")
  } 
  
  tmp <- data[vars]
  tmp <- tmp[do.call(order, tmp), ]
  lag <- ntimes2span - 1
  xTime <- convertTimes2numeric(tmp[[times]])
  miss.days <- sort(unique(xTime))[1:lag]
  
  #time.diffs are available - check them
  if (avail.times.diffs && 
      !all(is.na(data[[times.diffs]][xTime %in% miss.days])))
    stop("The times.diffs column in data does not have the appropriate missing values for the ", 
         "inital times of each individual\n",
         "Set avail.times.diffs to FALSE to have 'byIndv4Intvl_GRsDiff' calculate them")
  
  if (any(is.na(data[[times]])))
    warning(paste("Some values of ",times,
                  " are missing, which can result in merge producing a large data.frame", 
                  sep = ""))
  if (any(unlist(lapply(as.list(data[individuals]), 
                        function(f)
                          any(is.na(f))))))
    warning(paste("Some values of the factors in individuals are missing, ",
                  "which can result in merge producing a large data.frame", sep = ""))
  
  #Form time differences in a way that every first time is the same Time
  # - setting first time point to missing results in the growth rates also being NA
  if (!avail.times.diffs)
  { 
    #nspan (t)  2   3    4  5    6  7
    #nmiss      1   2    3  4    5  6 = lag
    #posn1      2   2    3  3    4  4 = ceiling((t+1)/2)
    #(t+1)/2  1.5   2  2.5  3  3.5  4
    #NA move    0   1    1  2    2  3 = nmiss - (posn1 -1) = floor((t+1) /2) - 1
    
#    tmp[times.diffs] <- calcLagged(xTime, operation ="-", lag = lag)
    tmp[times.diffs] <- xTime - unlist(by(xTime, tmp[individuals], FUN = calcLagged))
    
    #Make sure that first times.diffs are NA
    if (!all(is.na(tmp[[times.diffs]][xTime %in% miss.days])))
      tmp[[times.diffs]][xTime %in% miss.days] <- NA
    rownames(tmp) <- NULL
  }
  
  #Calculate the water traits
  names(water.trait.names) <- which.trait.types
  #Form WU (because time.diffs is NA for first time, so will the growth rates)
  tmp.indv <- split(tmp, f = as.list(tmp[individuals]))
  WU <- lapply(tmp.indv, function(d)
  {
    if (is.null(water.added))
      WU <- WU(weight.after = d[[weight.after]], weight.before = d[[weight.before]],
               time.diffs = d[[times.diffs]], lag = lag)
    else
      WU <- WU(weight.after = d[[weight.after]], water.added = d[[water.added]],
               time.diffs = d[[times.diffs]], lag = lag)
    return(WU)
  })
  WU <- unsplit(WU, tmp[[individuals]])
  WU[tmp[[times.diffs]][xTime %in% miss.days]] <- NA
  
  if ("WU" %in% which.trait.types)
    tmp[water.trait.names["WU"]] <- WU
  
  
  #Form WUR (because time.diffs is NA for first time, so will the growth rates)
  if (any(c("WUR", "WUI") %in% which.trait.types))
    WUR <- WU/tmp[[times.diffs]]
  if ("WUR" %in% which.trait.types)
    tmp[water.trait.names["WUR"]] <- WUR
  
  #Form WUI (because time.diffs is NA for first time, so will the growth rates)
  if ("WUI" %in% which.trait.types)
  { 
    if (is.null(responseAGR))
      stop("WUI is in which.trait.types, but responseAGR is NULL")
    tmp[water.trait.names["WUI"]] <- WUI(tmp[[responseAGR]], WUR)
  }
  
  #Reposition the GRs, if necessary
  nmove <- floor((ntimes2span+1) /2) - 1
  if (nmove > 0)
  {
    if (!avail.times.diffs)
      cols2move <- c(times.diffs, water.trait.names)
    else
      cols2move <- water.trait.names
    tmp <- split(tmp, f = as.list(tmp[individuals]))
    tmp <- lapply(tmp,
                  function(x, cols2move, nmove)
                  {
                    x[cols2move] <- x[c((nmove+1):nrow(x), 1:nmove), cols2move]
                    return(x)
                  },
                  cols2move = cols2move, nmove = nmove)
    tmp <- do.call(rbind, tmp)
  }
  
  #Remove NAs in individuals and time.factor in tmp
  if (any(is.na(tmp[[times]])))
    tmp <- tmp[!is.na(tmp[[times]]), ]
  if (any(unlist(lapply(as.list(tmp[individuals]), 
                        function(f)
                          any(is.na(f))))))
    for (f in individuals)
      tmp <- tmp[!is.na(tmp[f]),]
  #Keep times.diffs in data if they were used
  if (avail.times.diffs)
    tmp <- tmp[,-match(times.diffs, names(tmp))]
  #Remove unused times.diffs from data so time.diffs used are added from tmp
  if (!avail.times.diffs && times.diffs %in% names(data))
    data <- data[,-match(times.diffs, names(data))]
  if (any(water.trait.names %in% names(data)))
    data <- data[,-na.omit(match(water.trait.names, names(data)))]
  data <- dplyr::left_join(data, tmp, by = intersect(names(tmp), names(data)))
  data  <- data[do.call(order, data),]
  return(data)
}

##### I need to check that get missing values at the start of each Snapshot.ID.Tag.
"byIndv4Times_periodicRates" <- function(data, responses = NULL, 
                                         individuals = "Snapshot.ID.Tag", 
                                         times = "DAP", 
                                         columns2duplicate = NULL, 
                                         reqd.times.diff = 1, 
                                         avail.times.diffs = FALSE)
{ 
  #Check that responses, individuals, times and, if avail.times.diffs == TRUE, tims.diffs are in data
  times.diffs <- paste(times, "diffs", sep=".")
  vars <- c(individuals, times, columns2duplicate, responses)
  if (avail.times.diffs) 
    vars <- c(vars, times.diffs)
  checkNamesInData(vars, data = data)
  
  tmp <- data[vars]
  tmp <- tmp[do.call(order, tmp), ]
  lag <- 1
  xTim <- convertTimes2numeric(tmp[[times]])
  miss.days <- sort(unique(xTim))[1:lag]
  
  #Check for missing values that can mess up merging
  if (any(is.na(data[[times]])))
    warning(paste("Some values of ",times,
                  " are missing, which can result in merge producing a large data.frame", 
                  sep = ""))
  if (any(unlist(lapply(as.list(data[individuals]), 
                        function(f)
                          any(is.na(f))))))
    warning(paste("Some values of the factors in individuals are missing, ",
                  "which can result in merge producing a large data.frame", sep = ""))
  
  #time.diffs are available - check them
  if (avail.times.diffs) 
  { 
    if (!all(is.na(data[[times.diffs]][xTim %in% miss.days])))
      stop("The times.diffs column in data does not have the appropriate missing values for the ", 
           "inital times of each individual\n",
           "Set avail.times.diffs to FALSE to have 'byIndv4Intvl_GRsDiff' calculate them")
    tmp <- cbind(data, xTim)
  } else
  { 
  #Form time differences in a way that every first time is the same Time
  # - setting first time point to missing results in the growth rates also being NA
    #nspan (t)  2   3    4  5    6  7
    #nmiss      1   2    3  4    5  6 = lag
    #posn1      2   2    3  3    4  4 = ceiling((t+1)/2)
    #(t+1)/2  1.5   2  2.5  3  3.5  4
    #NA move    0   1    1  2    2  3 = nmiss - (posn1 -1) = floor((t+1) /2) - 1
    d <- cbind(tmp, xTim)
    indv <- as.list(d[individuals])
    d <- split(d, indv)
    d <- lapply(d, 
                function(df)
                {  
                  df[times.diffs] <- calcLagged(df[["xTim"]], operation ="-", lag = lag)
                  return(df)
                })
    tmp <- unsplit(d, indv)
    
    #Make sure that first times.diffs are NA
    if (!all(is.na(tmp[[times.diffs]][xTim %in% miss.days])))
      tmp[[times.diffs]][xTim %in% miss.days] <- NA
    rownames(tmp) <- NULL
  }
  if (!all.equal(tmp[times.diffs]/reqd.times.diff, round(tmp[times.diffs]/reqd.times.diff)))
    stop("Some of the values in the times.diffs column in data are not integer multiples of reqd.times.diff")

  #Check whether there all the time.diffs are already equal to reqd.times.diff required
  if (all(tmp[times.diffs][!is.na(tmp[times.diffs])] == reqd.times.diff, na.rm = TRUE))
  { 
    warning("All time differences are equal to reqd.times.diff for ", times, "; no action taken.")
    full.dat <- data
  } else #have time diffs not equal to rates.times
  {
    min.time <- min(xTim, na.rm = TRUE)
    max.time <- max(xTim, na.rm = TRUE)
    xtimes <- seq(min.time, max.time, by = reqd.times.diff)
    indv <- levels(factor(tmp[[individuals]]))
    full.dat <- data.frame()
    full.dat <- data.frame(xTim = rep(xtimes, times = length(indv)))
    full.dat <- cbind(dae::fac.gen(list(individ = indv), each = length(xtimes)), full.dat)
    if (is.factor(tmp[[times]]))
    { 
      full.dat[times] <- factor(full.dat$xTim) #add factor fpr times
      names(full.dat)[c(1,3)] <- c(individuals, times)
      xtimes <- "xTim"
      full.dat <- dplyr::left_join(full.dat, tmp, by = c(individuals, times, xtimes))
    } else
    {
      names(full.dat) <- c(individuals, times)
      xtimes <- times
      full.dat <- dplyr::left_join(full.dat, tmp, by = c(individuals, times))
    }
    #Now propagate Rates for DAP.diffs neq reqd.times.diff - needs to be done independently for each individual 
    indv.full <- as.list(full.dat[individuals])
    full1 <- split(full.dat, indv.full)
    propagateRates <- function(df, t, responses, columns2duplicate) 
    {
      time <- df[[xtimes]][t]
      time.diff <- df[[times.diffs]][t]
      to <- seq((time-(time.diff-1)), time-1, by = reqd.times.diff)
      for (kresp in responses)
      {
        df[[kresp]][df[[times]] %in% to] <- df[[kresp]][df[[times]] == time]
      # This makes the times.diffs applicable to the WUR, not the WU that remains.
      # However, the times.diffs as they stand are more apt for the WU than the equally spaced time 
      #   differences and so should be kept.
#        df[[times.diffs]][df[[times]] %in% c(time,to)] <- reqd.times.diff
      }
      if (!is.null(columns2duplicate))
      { 
        for (kcol in columns2duplicate)
          df[[kcol]][df[[times]] %in% to] <- df[[kcol]][df[[times]] == time]
      }
      return(df)
    }
    
    f1 <- mapply(function(full)
    {
      if (any(full[times.diffs] > reqd.times.diff))
      {
        which.bigger <- which(full[times.diffs] > reqd.times.diff )
        for (t in which.bigger)
          full <- propagateRates(full, t, responses, columns2duplicate) 
      }
      return(full)
    }, full = full1, SIMPLIFY = FALSE)
    full.dat <- unsplit(f1, indv.full)
    full.dat <- full.dat[vars]
    data <- data[, -match(intersect(c(times.diffs, columns2duplicate, responses), 
                                    names(data)), 
                          names(data))]
    full.dat <- dplyr::left_join(full.dat[vars], data, by = c(individuals, times))
    full.dat <- full.dat[c(setdiff(names(full.dat), responses), responses)]
  }
  return(full.dat)
}


