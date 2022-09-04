#Various functions for doing calculations.

#Function to calculate the cumulative sum, ignoring the first element if exclude.1st is TRUE
"cumulate" <- function(x, exclude.1st=FALSE, na.rm = FALSE, ...)
{ 
  if (na.rm)
    x <- x[!(is.na(x))]
  sum <- x
  if (exclude.1st)
    sum[-1] <- cumsum(x[-1])
  else
    sum <- cumsum(x)
  return(sum)
}

#Functions to calculate growth rates between successive imagings
"AGRdiff" <- function(x, time.diffs, lag=1)
{ 
  x.diffs <- calcLagged(x, operation = "-", lag = lag)
  x.diffs <- x.diffs / time.diffs
}
"PGR" <- function(x, time.diffs, lag=1)
{ 
  x.rates <- calcLagged(x, operation = "/", lag = lag)
  x.rates <- x.rates ^ (1/time.diffs)
}

"RGRdiff" <- function(x, time.diffs, lag=1)
{ 
  x.rates <- log(PGR(x, time.diffs, lag = lag))
}

"WUI" <- function(response, water)
{ 
  response.WUI <- ifelse(water != 0, 
                         response / water, 
                         NA)
  return(response.WUI)
}

#Functions to calculate a single-valued function, including the observation has the value of the function

"which.funct.value" <- function(x, FUNCT = NULL, na.rm = TRUE, ...)
{ 
  funct <- get(FUNCT)
  funct <- match.fun(funct)
  #  k <- which(x==funct(x, ...))
  if (all(is.na(x)))
    k <- NA
  else
  { 
    adiff <- abs(x - funct(x, na.rm = na.rm, ...))
    k <- which(adiff == min(adiff, na.rm = na.rm))
    if (length(k) == 0)
      k <- NA
    else if (length(k) > 1)
      k <- k[1]
  }
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

#Function to get the values for specified times and, optionally, a column with the times factor 
"getTimesSubset" <- function(data, responses, 
                             individuals = "Snapshot.ID.Tag", times = "DAP", 
                             which.times, suffix = NULL, sep.suffix.times = ".", 
                             include.times = FALSE, include.individuals = FALSE)
{ 
  n <- dim(data)[1]
  subset <- rep(FALSE, n)
  for (day in which.times)
  { 
    this.times <- data[times] == day
    subset <- subset | this.times
  }
  data.sub <- data[subset, ]
  data.sub <- data.sub[do.call(order, data.sub), ]
  if (!is.null(suffix))
    new.responses <- unlist(lapply(responses, 
                                   function(name, suffix){paste(name, suffix, sep=sep.suffix.times)}, 
                                   suffix=suffix))
  else
    new.responses <- responses
  if (include.individuals)
  {
    responses <- c(responses,individuals)
    new.responses <- c(new.responses,individuals)
  }
  if (include.times)
  {
    responses <- c(responses,times)
    if (!is.null(suffix))
      new.responses <- c(new.responses, paste(times, suffix, sep = sep.suffix.times))
    else
      new.responses <- c(new.responses,times)
  }
  data.sub <- data.sub[responses]
  names(data.sub) <- new.responses 
  return(data.sub)
}

#Functions to do calculations between successive dates 
# - does not assume same number time points for all individuals
#"Replace"  <- function(x, y) {z <- y}
"calcLagged" <- function(x, operation = NULL, lag=1)
  #This function replaces the observations with values calculated  
  # (i) for positive lag, itself and the value lag observations before it, 
  # (ii) for negative lag, itself and the value lag observations after it.
  #operation specifies calculation to be made on the pair of  values 
  #It returns as many values as are in data, the 1st lag values being NA
{ 
  n <- length(x)
  nl <- n-abs(lag)
  if (is.null(operation))
  { 
    if (lag > 0)
      x[(lag+1):n] <- x[1:nl]
    else
    { 
      if (lag < 0)
        x[1:nl] <- x[(abs(lag)+1):n]
    }
  }
  else
  { 
    FUN <- get(operation)
    FUN <- match.fun(FUN)
    if (lag > 0)
      x[(lag+1):n] <- FUN(x[(lag+1):n], x[1:nl])
    else
    { 
      if (lag < 0)
        x[1:nl] <- FUN(x[1:nl], x[(abs(lag)+1):n])
      else
        x[1:n] <- FUN(x[1:n], x[1:n])
    }
  }
  if (lag > 0)
    x[1:lag] <- NA
  else
    x[(nl+1):n] <- NA
  return(x)
}

#Function to test if any values in a set of values are anomalous
# in being outside specified limits 
"anom" <- function(x, lower=NULL, upper=NULL, na.rm = TRUE)
{ if (is.null(lower))
{ if (is.null(upper))
  stop("Must supply at least a lower or an  upper limit below or above which values are anomalous")
  else
    anom <- any(x > upper, na.rm=na.rm)
} else
{ if (is.null(upper))
  anom <- any(x < lower, na.rm=na.rm)
else
  anom <- any(x < lower, na.rm=na.rm) || any(x > upper, na.rm=na.rm)
}
  return(anom)
}

"twoLevelOpcreate" <- function(data, responses, treatment.factor="Treatment.1", 
                               suffices.treatment=c("Cont","Salt"), control=1, 
                               columns.suffixed = NULL, 
                               operations = "/", suffices.results="OST", 
                               columns.retained = c("Snapshot.ID.Tag",
                                                    "Smarthouse","Lane", 
                                                    "Zone","cZone","SHZone",
                                                    "ZLane","ZMainunit","cMainPosn", 
                                                    "Genotype.ID"),
                               by = c("Smarthouse","Zone","ZMainunit"))
{ 
  #Check that response and treatment.factor are in data
  checkNamesInData(c(responses, treatment.factor), data)
  
  trt.lev <- levels(data[[treatment.factor]])
  ntrt <- length(trt.lev)
  if (ntrt != 2)
    stop("The treatment factor must have only two levels")
  if (!(control %in% c(1,2)))
    stop("contol must be either 1 or 2")
  if (!all(by %in% columns.retained))
    stop("The names in by must be included in columns.retained")
  if (treatment.factor %in% columns.retained)
    stop("The treatment.factor should not be included in columns.retained")
  treated <- c(1,2)[c(1,2)!=control]
  if (!is.null(columns.suffixed))
    if (any(!(unique(columns.suffixed) %in% c(responses,columns.retained))))
      stop("The columns.suffixed must be included in columns.retained")
  
  #identify the functions specified in operations and check the number of them
  nresp <-  length(responses)
  fun <- vector(mode="list", length=nresp)
  fun <- lapply(operations,
                function(operations)
                  fun <- match.fun(get(operations))
  )
  if (length(operations)==1 && nresp != 1)
    fun[2:nresp] <- fun[1]
  if (length(fun) != nresp)
    stop("The length of operations is neither 1 nor the length of responses")
  if (length(suffices.results) == 1 && nresp != 1)
    suffices.results <- rep(suffices.results, nresp)
  if (length(suffices.results) != nresp)
    stop("The length of suffices.results is neither 1 nor the length of responses")
  
  #Split data according to the two treatments 
  uresponses <- unique(responses)
  ucols.suffixed <- unique(columns.suffixed)
  data <- data[c(treatment.factor,columns.retained, uresponses)]
  data <- split(data, data[[treatment.factor]])
  data[[1]] <- data[[1]][-match(treatment.factor,names(data[[1]]))]
  data[[2]] <- data[[2]][-match(treatment.factor,names(data[[2]]))]
  
  #Create data frame with the retained columns and columns for the GR values for each treatment
  data.out <- data[[treated]]
  names(data.out)[match(c(uresponses,ucols.suffixed), names(data.out))] <- 
    paste(c(uresponses,ucols.suffixed),suffices.treatment[treated],sep=".")
  data.out <- merge(data.out, data[[control]][c(by,uresponses,ucols.suffixed)], by=by,sort=FALSE)
  names(data.out)[match(c(uresponses,ucols.suffixed), names(data.out))] <- 
    paste(c(uresponses,ucols.suffixed),suffices.treatment[control],sep=".")
  
  #Add the results of the operations
  data.out[paste(responses,suffices.results,sep=".")] <- 
    lapply(1:nresp, 
           function(k, responses, fun, data)
           { out <- fun[[k]](data[[paste(responses[[k]],suffices.treatment[treated],sep=".")]], 
                             data[[paste(responses[[k]],suffices.treatment[control],sep=".")]])
           return(out)
           }, responses=responses, fun=fun, data=data.out)
  return(data.out)
}

