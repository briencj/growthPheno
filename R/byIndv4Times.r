#Functions to produce values for all times for eacfh individual
"byIndv4Times_GRsDiff" <- function(data, responses, 
                                   individuals = "Snapshot.ID.Tag", times = "DAP", 
                                   which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
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
      responses.GR <- paste(responses, "AGR", sep=".")
    else
      responses.GR <- paste(responses, suffices.rates[match("AGR",opt)], sep=".")
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                              FUN = AGRdiff, 
                                              time.diffs = tmp[[times.diffs]], lag = lag))
    responses.GRs <- c(responses.GRs, responses.GR)
  }
  
  #Form PGR (because time.diffs is NA for first time, so will the growth rates)
  if ("PGR" %in% opt)
  { 
    if (is.null(suffices.rates))
      responses.GR <- paste(responses, "PGR", sep=".")
    else
      responses.GR <- paste(responses, suffices.rates[match("PGR",opt)], sep=".")
    tmp[responses.GR] <- as.data.frame(lapply(tmp[responses], 
                                              FUN = PGR, 
                                              time.diffs = tmp[[times.diffs]], lag = lag))
    responses.GRs <- c(responses.GRs, responses.GR)
  }
  
  #Form RGR (because time.diffs is NA for first time, so will the growth rates)
  if ("RGR" %in% opt)
  { 
    if (is.null(suffices.rates))
      responses.GR <- paste(responses, "RGR", sep=".")
    else
      responses.GR <- paste(responses, suffices.rates[match("RGR",opt)], sep=".")
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
    tmp <- split(tmp, f = as.list(tmp[individuals], simplify=FALSE))
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
  data <- merge(data, tmp, by = c(individuals, times), sort = FALSE, all.x = TRUE)
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
                          na.x.action, na.y.action, sep, ...)
{
  #Split data frame by each combination of the individuals factors
  old.names <- names(data)
  tmp <- split(data, as.list(data[individuals]), sep=sep)
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
  indices <- strsplit(indices, split=sep, fixed=TRUE)
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
                                      which.rates = c("AGR","RGR"), suffices.rates = NULL, 
                                      avail.times.diffs = FALSE, ntimes2span = 2, 
                                      extra.derivs = NULL, suffices.extra.derivs=NULL, 
                                      sep=".", 
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
                        rates = splrates, suffices.rates = splsuffices,
                        extra.derivs = extra.derivs, 
                        suffices.extra.derivs = suffices.extra.derivs, 
#                        deriv = derivs, suffices.deriv = suffices.derivs, 
#                        extra.rate = extra.rate, 
                        na.x.action = na.x.action, na.y.action = na.y.action, 
                        sep = sep, ...)
    if (rates.method == "differences")
    {
      smth <- byIndv4Times_GRsDiff(data = smth, response.smoothed, 
                                   individuals=individuals,times=times, 
                                   which.rates = grates, ntimes2span = ntimes2span, 
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
                           na.x.action=na.x.action, na.y.action=na.y.action, sep = sep, ...)
      if (rates.method == "differences" && ntimes2span != 2)
        smth <- byIndv4Times_GRsDiff(data = smth, response.smoothed, 
                                     individuals=individuals,times=times, 
                                     which.rates = grates, ntimes2span = ntimes2span, 
                                     avail.times.diffs = avail.times.diffs)
      subdat[times] <- convertTimesExnumeric(subdat[[times]], data[[times]])
      smth <- rbind(smth, subdat)
    }
    smth <- smth[do.call(order, smth), ]
    if (rates.method == "differences" && ntimes2span == 2)
      smth <- byIndv4Times_GRsDiff(data = smth, response.smoothed, 
                                   individuals=individuals,times=times, 
                                   which.rates = grates, ntimes2span = ntimes2span, 
                                   avail.times.diffs = avail.times.diffs)
  }
  
  return(smth)
}
