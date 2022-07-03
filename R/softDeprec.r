"splitContGRdiff" <- function(data, responses, 
                              individuals = "Snapshot.ID.Tag", INDICES = NULL, 
                              which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                              times.factor = "Days", avail.times.diffs = FALSE, 
                              ntimes2span = 2)
{ 
  .Deprecated(old = "splitContGRdiff", new = "byIndv4Times_GRsDiff", package = "growthPheno",
              msg = paste0("'splitContGRdiff' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4Times_GRsDiff' instead."))
  
  if (!is.null(INDICES))
    individuals <- INDICES
  
  options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  vars <- c(individuals, times.factor, responses)
  times.diffs <- paste(times.factor, "diffs", sep=".")
  if (avail.times.diffs)
  { 
    if (times.diffs %in% names(data))
      vars <- c(vars, times.diffs)
    else
      stop("The column ", times.diffs, ", expected to contain the times differences, is not in data")
  }
  
  #Get columns needed for calculation and order for individuals, then times.factor
  if (any(is.na(match(vars, names(data)))))
    stop("One or more of response, individuals and times.factor are not in data")
  tmp <- data[vars]
  tmp <- tmp[do.call(order, tmp), ]
  lag <- ntimes2span - 1
  xTime <- dae::as.numfac(tmp[[times.factor]])
  miss.days <- sort(unique(xTime))[1:lag]
  
  #time.diffs are available - check them
  if (avail.times.diffs && 
      !all(is.na(data[[times.diffs]][xTime %in% miss.days])))
    stop("The times.diffs column in data does not have the appropriate missing values for the ", 
         "inital days of each individual\n",
         "Set avail.times.diffs to FALSE to have 'splitContGRdiff' calculate them")
  
  if (any(is.na(data[[times.factor]])))
    warning(paste("Some values of ",times.factor,
                  " are missing, which can result in merge producing a large data.frame", 
                  sep = ""))
  if (any(unlist(lapply(as.list(data[individuals]), 
                        function(f)
                          any(is.na(f))))))
    warning(paste("Some values of the factors in individuals are missing, ",
                  "which can result in merge producing a large data.frame", sep = ""))
  
  #Form time differences in a way that every first day is the same Day
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
  #Form AGR (because Day.diffs is NA for first day, so will the growth rates)
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
  
  #Form PGR (because Day.diffs is NA for first day, so will the growth rates)
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
  
  #Form RGR (because Day.diffs is NA for first day, so will the growth rates)
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
                  function(tmp, cols2move, nmove)
                  {
                    tmp[cols2move] <- tmp[c((nmove+1):nrow(tmp), 1:nmove), cols2move]
                    return(tmp)
                  },
                  cols2move = cols2move, nmove = nmove)
    tmp <- do.call(rbind, tmp)
  }
  
  #Remove NAs in individuals and time.factor in tmp
  if (any(is.na(tmp[[times.factor]])))
    tmp <- tmp[!is.na(tmp[[times.factor]]), ]
  if (any(unlist(lapply(as.list(tmp[individuals]), 
                        function(f)
                          any(is.na(f))))))
    for (f in individuals)
      tmp <- tmp[!is.na(tmp[f]),]
  tmp <- tmp[,-match(responses,names(tmp))]
  if (avail.times.diffs)
    tmp <- tmp[,-match(times.diffs, names(tmp))]
  if (!avail.times.diffs && times.diffs %in% names(data))
    data <- data[,-match(times.diffs, names(data))]
  if (any(responses.GR %in% names(data)))
    data <- data[,-match(responses.GRs, names(data))]
  data <- merge(data, tmp, by = c(individuals, times.factor), sort = FALSE, all.x = TRUE)
  data  <- data[do.call(order, data),]
  return(data)
}

#### The only function that uses fitSpline is splitSplines and so, when splitSplines is 
#### deprecated, fitSpline should also be deprecated.
#Function to fit a spline using smooth.spline or JOPS
"fitSpline" <- function(data, response, response.smoothed, x, 
                        smoothing.method = "direct", 
                        spline.type = "NCSS",  df=NULL, lambda = NULL, 
                        npspline.segments = NULL, correctBoundaries = FALSE, 
                        deriv = NULL, suffices.deriv = NULL, extra.rate = NULL, 
                        na.x.action = "exclude", na.y.action = "trimx", ...)
{ 
  #This table has been made obsolete with the introduction of extra.rate
  #  Result   deriv   suffix.deriv   AGR     RGR
  # direct smoothing
  # AGR, RGR    1          AGR      NULL     RGR    
  #   AGR       1          AGR      NULL    NULL
  #   RGR       1         NULL      NULL     RGR
  # log-smoothing
  # AGR, RGR    1          RGR       AGR    NULL
  #   AGR       1         NULL       AGR    NULL
  #   RGR       1          RGR      NULL    NULL
  
  if (!all(c(response,x) %in% names(data)))
    stop(paste("One or more of", response, "and", x, "is missing from ", 
               deparse(substitute(data))))
  
  #check input arguments
  impArgs <- match.call()
  if ("na.rm" %in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  if ("smoothing.scale" %in% names(impArgs))
    stop("smoothing.scale has been deprecated; use smoothing.method")
  
  #Check that required cols are in data
  checkNamesInData(c(response, x), data = data)
  
  if (missing(response.smoothed))
    stop("The argument response.smoothed has not been supplied.")
  
  smethods <- c("direct", "logarithmic")
  smethod <- smethods[check.arg.values(smoothing.method, options=smethods)]
  stype <- c("NCSS", "PS")
  stype <- stype[check.arg.values(spline.type, options=stype)]
  if (is.null(lambda) && stype == "PS")
    stop("Must specify lambda for spline.type set to PS")
  if (!is.null(df) && !is.null(lambda) && stype == "NCSS")
    stop("Only one of df and lambda can be specified for spline.type NCSS")
  
  na.x <- na.y <- c("exclude", "omit", "fail")
  na.y <- c(na.y, "allx", "trimx", "ltrimx", "utrimx")
  na.act.x <- na.x[check.arg.values(na.x.action, options=na.x)]
  na.act.y <- na.y[check.arg.values(na.y.action, options=na.y)]  
  
  if (correctBoundaries & !is.allnull(deriv))
    stop("Unable to estimate derivatives when correctBoundaries = TRUE")
  
  if (!is.null(deriv) & !is.null(suffices.deriv))
    if (length(deriv) != length(suffices.deriv))
      stop("The number of names supplied must equal the number of derivatives specified")
  if (!is.null(extra.rate)) 
  {
    if (!(1 %in% deriv))
      stop("To form an extra rate, 1 must be included in deriv so that the first derivative is obtained")
    else
      kagr <- match(1, deriv)
    if (smethod == "direct" && extra.rate != "RGR")
      stop("Only the RGR can be formed from the first derivative when direct smoothing is used")
    if (smethod == "logarithmic" && extra.rate != "AGR")
      stop("Only the AGR can be formed from the first derivative when logarithmic smoothing is used")
    if (is.null(names(extra.rate)))
      names(extra.rate) <- extra.rate
  }
  
  tmp <- as.data.frame(data)
  #Transform data, if required
  if (smethod == "logarithmic")
    tmp[[response]] <- log(tmp[[response]])
  
  #Convert any infinite values to missing
  if (any(is.infinite(tmp[[response]])))
  {
    tmp[[response]][is.infinite(tmp[[response]])] <- NA
    warning("Some infinite values have been converted to missing")
  }
  
  #Set up fit data.frame
  if (is.null(response.smoothed))
    response.smoothed <- paste(response,"smooth",sep=".")
  fit.names <- c(x, response.smoothed)
  if (!is.allnull(deriv))
  {
    if (is.allnull(suffices.deriv))
      fit.names <- c(fit.names, paste(response.smoothed,".dv",deriv,sep=""))
    else
      fit.names <- c(fit.names, paste(response.smoothed, suffices.deriv, sep="."))
  }
  
  #Process missing values
  #tmp will have no missing values and is what is used in the fitting
  #data retains missing values and is used to obtain the returned data.frame
  #x.pred has the x-values for which predictions are required 
  #  - it should include all the x values that are returned by smooth.spline;
  #    it will not include any x values that are missing, but may include
  #    x values for which there are missing y values, depending on the settings 
  #    of na.y.action, and these x values will not have been supplied to smooth.spline.
  nobs <- nrow(tmp)
  if (nobs == 0)
  {
    warning("A response with no data values supplied")
  } else
  {
    if (na.act.x == "fail")
      stop("na.x.action is to set to fail and there are missing x values")
    else #remove any observations with missing x values
    {
      if (na.act.x %in% c("omit", "exclude"))
      {
        tmp <- tmp[!is.na(tmp[[x]]), ]
        if (na.act.x == "omit")
          data <- tmp
      }
    }
    x.pred <- tmp[[x]]
    if (na.act.y == "fail")
      stop("na.y.action is to set to fail and there are missing y values")
    else #Are there any missing response values now
    {
      if (na.act.y %in% c("omit", "exclude"))
      {
        x.pred <- tmp[!is.na(tmp[[response]]), ][[x]]
      } else
      {
        if (grepl("trimx", na.act.y, fixed = TRUE))
        {
          tmp <- tmp[order(tmp[[x]]), ]
          y.nonmiss <- c(1:nrow(tmp))[!is.na(tmp[[response]])]
          x.pred <- tmp[[x]]
          if (length(y.nonmiss) > 0)
          {
            if (na.act.y %in% c("trimx", "ltrimx"))
            {
              x.pred <- x.pred[y.nonmiss[1]:length(tmp[[x]])]
              y.nonmiss <- y.nonmiss - y.nonmiss[1] + 1
            }
            if (na.act.y %in% c("trimx", "utrimx"))
              x.pred <- x.pred[1:y.nonmiss[length(y.nonmiss)]]
            y.nonmiss <- diff(y.nonmiss)
            if (any(y.nonmiss >= 3))
              warning("There are runs of 3 or more contiguous y values that are missing")
          } else
            x.pred <- NULL
        }
      }
      tmp <- tmp[!is.na(tmp[[response]]), ]
      if (na.act.y == "omit")
        data <- tmp
    }
  }
  
  #What are the distinct x values
  distinct.xvals <- sort(unique(tmp[[x]]))
  tol <- 1e-06 * IQR(distinct.xvals)
  distinct.xvals <- remove.repeats(distinct.xvals, tolerance = tol)
  if (length(distinct.xvals) < 4) 
  { 
    #< 4 distinct values and so all fitted values are set to NA
    warning(paste("Need at least 4 distinct x values to fit a spline",
                  "- all fitted values set to NA", sep = " "))
    #Set up fit data.frame
    if (is.null(response.smoothed))
      response.smoothed <- paste(response,"smooth",sep=".")
    fit.names <- c(x, response.smoothed)
    if (!is.allnull(deriv))
    {
      if (is.allnull(suffices.deriv))
        fit.names <- c(fit.names, paste(response.smoothed,".dv",deriv,sep=""))
      else
        fit.names <- c(fit.names, paste(response.smoothed, suffices.deriv, sep="."))
    }
    #Add extra,rate if required
    if (!is.null(extra.rate))
      fit.names <- c(fit.names, paste(response.smoothed, names(extra.rate), sep="."))
    fit <- as.data.frame(matrix(NA, nrow=nrow(data), ncol = length(fit.names)))
    colnames(fit) <- fit.names
    fit[x] <- data[[x]]
    fit.spline <- NULL
  } else
  { 
    #smooth and obtain predictions corresponding to x.pred
    fitcorrectBoundaries <- correctBoundaries
    if (stype == "NCSS" && length(distinct.xvals) <= 5 && correctBoundaries)
    {
      warning(paste("Need more than 5 distinct x values to correct the end-points of a spline",
                    "- no corrections made", sep = " "))
      fitcorrectBoundaries <- FALSE
    }
    if (stype == "NCSS")
    {
      if (is.null(df))
      {
        if (is.null(lambda))
          fit.spline <- ncsSpline(tmp[c(x, response)], correctBoundaries = fitcorrectBoundaries, 
                                  ...)
        else
          fit.spline <- ncsSpline(tmp[c(x, response)], correctBoundaries = fitcorrectBoundaries, 
                                  lambda = lambda, ...)
      } else
        if (is.null(lambda))
          fit.spline <- ncsSpline(tmp[c(x, response)], correctBoundaries = fitcorrectBoundaries, 
                                  df = df, ...)
    } else #PS
    {
      #Determine npspline.segments for a full set of x
      if (is.null(npspline.segments))
        npspline.segments <- max(10, ceiling((nrow(data)-1)/2))
      #Adjust for the missing values in x
      nfit <- length(distinct.xvals)
      if (nfit < nobs)
      {
        nperseg <- nobs/npspline.segments
        npspline.segments <- ceiling(nfit/nperseg)
      }
      fit.spline <- pSpline(tmp[c(x, response)], npspline.segments = npspline.segments, lambda = lambda, ...)
    }
    x.pred <- remove.repeats(sort(x.pred), tolerance = tol)
    fit <- NULL
    if (length(x.pred) == length(fit.spline$x))
    {
      if (all(abs(x.pred - fit.spline$x) < tol))
      {  
        if (stype == "NCSS")
          fit <- list(fit.spline$x, fit.spline$y)
        else
          fit <- list(fit.spline$x, fit.spline$y)
      }
    }
    #Need to refit for current x.pred
    if (is.null(fit))
    {
      if (stype == "NCSS")
        fit <- predict.ncsSpline(fit.spline, x = x.pred, 
                                 correctBoundaries = fitcorrectBoundaries)
      else
        fit <- predict.pSpline(fit.spline, x = x.pred, npspline.segments = fit.spline$uncorrected.fit$npspline.segments)
    }
    rsmooth <- response.smoothed
    names(fit) <- c(x, rsmooth)
    #backtransform if transformed
    if (smethod == "logarithmic")
      fit[[rsmooth]] <- exp(fit[[rsmooth]])
    
    #get derivatives if required
    if (!correctBoundaries & !is.null(deriv))
    {
      for (d in deriv)
      {
        if (is.null(suffices.deriv))
          rsmooth.dv <- paste0(response.smoothed,".dv",d)
        else
        { 
          k <- match(d, deriv)
          rsmooth.dv <- paste(response.smoothed, suffices.deriv[k], sep=".")
        }
        if (stype == "NCSS")
          fit[[rsmooth.dv]] <- predict(fit.spline$uncorrected.fit, x = x.pred, deriv=d)$y
        else
          fit[[rsmooth.dv]] <- predict.pSpline(fit.spline, x = x.pred, 
                                               npspline.segments = fit.spline$uncorrected.fit$npspline.segments, 
                                               deriv=d)$y
      }
      #Add RGR if required
      if (!is.null(extra.rate) && extra.rate == "RGR")
      { 
        #Check have the required computed derivative 
        if (is.null(suffices.deriv))
          rsmooth.dv <- paste(response,".smooth.dv",1,sep="")
        else
        { 
          k <- match(1, deriv)
          rsmooth.dv <- paste(response.smoothed, suffices.deriv[k], sep=".")
        }
        if (!(rsmooth.dv %in% names(fit)))
          stop("First derivative not available to calculate RGR")
        fit[[paste(rsmooth,names(extra.rate),sep=".")]] <- fit[[rsmooth.dv]]/fit[[rsmooth]]
      }
      #Add AGR if required
      if (!is.null(extra.rate) && extra.rate == "AGR")
      { 
        #Check have the required computed derivative 
        if (is.null(suffices.deriv))
          rsmooth.dv <- paste(response,".smooth.dv",1,sep="")
        else
        { 
          k <- match(1, deriv)
          rsmooth.dv <- paste(response.smoothed, suffices.deriv[k], sep=".")
        }
        #get the extra derivative
        if (!(rsmooth.dv %in% names(fit)))
          stop("First derivative not available to calculate AGR")
        fit[[paste(rsmooth,names(extra.rate),sep=".")]] <- fit[[rsmooth.dv]]*fit[[rsmooth]]
      }
    }
    fit <- as.data.frame(fit)
    #Merge data and fit, preserving x-order in data
    x.ord <- order(data[[x]])
    fit <- merge(data[c(x,response)],fit, all.x = TRUE, sort = FALSE)
    fit <- fit[,-match(response, names(fit))]
    fit <- fit[order(fit[[x]]),]
    fit[x.ord,] <- fit
  }
  rownames(fit) <- NULL
  return(list(predictions = fit, fit.spline = fit.spline))    
}

#Fit splines to smooth the longitudinal trends in a set of individuals for a response
#Specify responses to be smoothed and then loop over the individuals
"indivSplines" <- function(data, 
                           response, response.smoothed, x, individuals, 
                           stype, smethod, df, npspline.segments, lambda, 
                           correctBoundaries, 
                           deriv, suffices.deriv, extra.rate, sep, 
                           na.x.action, na.y.action, ...)
{
  #Split data frame by each combination of the individuals factors
  old.names <- names(data)
  tmp <- split(data, as.list(data[individuals]), sep=sep)
  #Fit splines for each combination of the individuals factors
  tmp <- lapply(tmp, fitSpline, 
                response=response, response.smoothed=response.smoothed, 
                x = x, 
                spline.type = stype, smoothing.method = smethod, 
                df=df, npspline.segments = npspline.segments, lambda = lambda, 
                correctBoundaries = correctBoundaries, 
                deriv=deriv, suffices.deriv=suffices.deriv, extra.rate = extra.rate, 
                na.x.action=na.x.action, na.y.action=na.y.action, ...)
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
  tmp.smooth <- names(tmp)[-match(c(individuals,x), names(tmp))]
  tmp.smooth <- na.omit(match(tmp.smooth, names(data)))
  if (length(tmp.smooth) > 0)
    data <- data[ ,-tmp.smooth]
  tmp <- tmp[!is.na(tmp[[x]]), ]
  data <- merge(data, tmp, all.x = TRUE, sort=FALSE)
  #Rearrange columns so original column are followed by new columns
  new.names <- names(data)
  new.names <- new.names[-match(old.names, new.names)]
  data <- data[c(old.names, new.names)]
  rownames(data) <- NULL
  
  return(data)
}

#The function that controls the fitting
"splitSplines" <- function(data, response, response.smoothed = NULL, x,  
                           individuals = "Snapshot.ID.Tag", INDICES = NULL,
                           smoothing.method = "direct", smoothing.segments = NULL, 
                           spline.type = "NCSS", df=NULL, lambda = NULL, 
                           npspline.segments = NULL, 
                           correctBoundaries = FALSE, 
                           deriv = NULL, suffices.deriv=NULL, extra.rate = NULL, 
                           sep=".", 
                           na.x.action="exclude", na.y.action = "exclude", ...)
{ 
  .Deprecated(old = "splitSplines", new = "byIndv4Times_SplinesGRs", package = "growthPheno",
              msg = paste0("'splitSplines' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4Times_SplinesGRs' instead."))
  
  if (!is.null(INDICES))
    individuals <- INDICES
  
  if (!all(c(response,x) %in% names(data)))
    stop(paste("One or more of", response, "and", x, "is missing from ", 
               deparse(substitute(data))))
  
  impArgs <- match.call()
  if ("na.rm" %in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  if ("smoothing.scale" %in% names(impArgs))
    stop("smoothing.scale has been deprecated; use smoothing.method")
  
  smethods <- c("direct", "logarithmic")
  smethod <- smethods[check.arg.values(smoothing.method, options=smethods)]
  stype <- c("NCSS", "PS")
  stype <- stype[check.arg.values(spline.type, options=stype)]
  if (is.null(lambda) && stype == "PS")
    stop("Must specify lambda for spline.type set to PS")
  if (!is.null(df) && !is.null(lambda))
    stop("One of df and lambda must be NULL")
  
  if (is.null(response.smoothed))
    response.smoothed <- paste0(response, ".smooth")
  
  #Check extra.rate
  if (!is.null(extra.rate)) 
  {
    if (!(1 %in% deriv))
      stop("To form an extra rate, the first derivative must be obtained by incuding 1 in deriv")
    else
      kagr <- match(1, deriv)
    if (smethod == "direct" && extra.rate != "RGR")
      stop("Only the RGR can be formed from the first derivative when direct smoothing is used")
    if (smethod == "logarithmic" && extra.rate != "AGR")
      stop("Only the AGR can be formed from the first derivative when logarithmic smoothing is used")
  }
  
  #Check npspline.segments
  if (length(npspline.segments) > 1)
  { 
    if (is.null(smoothing.segments))
      stop("npspline.segments must be of length one in an unsegmented spline fit")
    else
    {
      if (length(npspline.segments) != length(smoothing.segments))
        stop("the number of values of npspline.segments should be one or ",
             "equal to the number of segments in a segmented spline fit")
    }
    if(!all(diff(unlist(smoothing.segments)) > 0))
      stop("the smoothing.segments are not a set of non-overlapping, successive intervals")
  }
  
  #Fit the splines to the individuals
  if (is.allnull(smoothing.segments))
    smth <- indivSplines(data= data,
                         response=response, response.smoothed=response.smoothed, 
                         x = x, individuals = individuals, 
                         stype = stype, smethod = smethod, 
                         df=df, npspline.segments = npspline.segments, lambda = lambda, 
                         correctBoundaries = correctBoundaries, 
                         deriv=deriv, suffices.deriv=suffices.deriv, extra.rate = extra.rate, 
                         sep = sep, 
                         na.x.action=na.x.action, na.y.action=na.y.action, ...)
  else    
  {
    knseg <- npspline.segments[1]
    smth <- data.frame()
    for (k in 1:length(smoothing.segments))
    {
      segm <- smoothing.segments[[k]]
      subdat <- data[(data[x] >= segm[1]) & (data[x] <= segm[2]),] 
      if (length(npspline.segments) > 1) knseg <- npspline.segments[k]
      smth <- rbind(smth, 
                    indivSplines(data = subdat, 
                                 response=response, response.smoothed=response.smoothed, 
                                 x = x, individuals = individuals, 
                                 stype = spline.type, smethod = smoothing.method, 
                                 df=df, npspline.segments = npspline.segments, lambda = lambda, 
                                 correctBoundaries = correctBoundaries, 
                                 deriv=deriv, suffices.deriv=suffices.deriv, extra.rate = extra.rate, 
                                 sep = sep, 
                                 na.x.action=na.x.action, na.y.action=na.y.action, ...))
    }
    smth <- smth[do.call(order, smth), ]
  }
  return(smth)
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
  .Deprecated(old = "splitValueCalculate", new = "byIndv4_ValueCalc", package = "growthPheno",
              msg = paste0("'splitValueCalculate' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4_ValueCalc' instead."))
  
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

#Functions for calculating derived responses
#Function to form the growth rates over an interval for a set of responses
"intervalGRdiff" <- function(responses, individuals = "Snapshot.ID.Tag", 
                             which.rates = c("AGR","PGR","RGR"), suffices.rates=NULL, 
                             times = "Days", start.time, end.time, 
                             suffix.interval, data)
{ 
  .Deprecated(old = "intervalGRdiff", new = "byIndv4Intvl_GRsDiff", package = "growthPheno",
              msg = paste0("'intervalGRdiff' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4Intvl_GRsDiff' instead."))
  
  options <- c("AGR","PGR","RGR")
  opt <- options[unlist(lapply(which.rates, check.arg.values, options=options))]
  if (!is.null(suffices.rates) & length(opt) != length(suffices.rates))
    stop("The length of of which.rates and suffices.rates should be equal")
  
  if (!all(sapply(list(start.time, end.time, suffix.interval), 
                  function(x) length(x) == 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length one")
  
  #Check that individuals and times are in data
  if (!all(c(individuals,times) %in% names(data)))
    stop("Indivduals and/or times are not in data")

  interval.resp <- cbind(getTimesSubset(responses, data = data, which.times = start.time, 
                                        times = times, 
                                        include.times = TRUE, 
                                        suffix = "start"),
                         getTimesSubset(responses, data = data, which.times = end.time, 
                                        times = times, 
                                        include.times = TRUE, 
                                        suffix = "end"))
  times.start <- paste(times,"start",sep=".")
  times.end <- paste(times,"end",sep=".")
  
  interval.resp[times.start] <- convertTimes2numeric(interval.resp[[times.start]])
  interval.resp[times.end] <- convertTimes2numeric(interval.resp[[times.end]])
  growth.rates <- lapply(responses, 
                         function(name, interval.resp, start.time, end.time, 
                                  which.rates, times = "Days")
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
                           { k <- k + 1
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
                           { k <- k +1
                             if ("RGR" %in% which.rates)
                               PGR <- exp(RGR)
                             else
                             { PGR <- (log(interval.resp[end.name]) - log(interval.resp[start.name])) / 
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
  growth.rates <- cbind(growth.rates, 
                        getTimesSubset(individuals, data = data, times = times, 
                                       which.times = start.time, suffix = NULL))
  return(growth.rates)
}

#Function to form the growth rates over an interval for a set of responses by averaging growth rates
"intervalGRaverage" <- function(responses, individuals = "Snapshot.ID.Tag", 
                                which.rates = c("AGR","RGR"), suffices.rates=c("AGR","RGR"), 
                                times = "Days", start.time, end.time, suffix.interval, 
                                data, sep=".", na.rm=TRUE)
{  
  .Deprecated(old = "intervalGRaverage", new = "byIndv4Intvl_GRsAvg", package = "growthPheno",
              msg = paste0("'intervalGRdiff' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4Intvl_GRsAvg' instead."))
  
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
                                       stop("Growth rates for at least some responses are not in data")
                                     return(response.grates)
                                   },
                                   responses=responses))
  #Check that individuals and times are in data
  if (!all(c(individuals,times) %in% names(data)))
    stop("Indivduals and/or times are not in data")
  #Get data for the times 
  times.vals <- unique(convertTimes2numeric(data[[times]]))
  times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
  interval.resp <- getTimesSubset(response.grates, data = data, which.times = times.vals, 
                                  times = times, include.times = TRUE)
  interval.resp <- cbind(getTimesSubset(individuals, data = data, which.times = times.vals, 
                                        times = times, 
                                        include.times = FALSE),
                         interval.resp)
  interval.resp[times] <- convertTimes2numeric(interval.resp[[times]])
  #calculate the weights
  interval.resp <- split(interval.resp, as.list(interval.resp[individuals]), sep=sep)
  interval.resp <- lapply(interval.resp, 
                          function(data, times = "Days")
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
                               rates <- splitValueCalculate(name.AGR, weights=times, 
                                                            individuals = individuals, 
                                                            FUN = "weighted.mean", data=interval.resp, 
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
                                 rates <- splitValueCalculate(name.RGR, weights=times, 
                                                              individuals = individuals, 
                                                              FUN = "weighted.mean", 
                                                              data=interval.resp, 
                                                              na.rm=na.rm, sep=sep)
                               else
                                 rates <- merge(rates, 
                                                splitValueCalculate(name.RGR, 
                                                                    weights=times, 
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
                                     times = "Days", start.time=NULL, end.time=NULL, 
                                     suffix.interval=NULL, data, sep=".", na.rm=TRUE, ...)
{  
  .Deprecated(old = "intervalValueCalculate", new = "byIndv4Intvl_ValueCalc", package = "growthPheno",
              msg = paste0("'intervalValueCalculate' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4Intvl_ValueCalc' instead."))
  
  #Trap which.levels and give message to replace it with which.values is not set
  tempcall <- list(...)
  if (length(tempcall) && "which.levels" %in% names(tempcall))
    stop("replace which.levels with which.values")
  
  if (!all(sapply(list(start.time, end.time, suffix.interval), 
                  function(x) length(x) <= 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length one or NULL")
  
  #Check that response is in data
  if (!all(c(response,individuals,times) %in% names(data)))
    stop("Some of the columns for response, indivduals and times are not in data")
  #Get data for the times
  if (all(is.null(c(start.time, end.time))))
    interval.resp <- data[c(individuals,response,times)]
  else
  { 
    times.vals <- unique(convertTimes2numeric(data[[times]]))
    if (is.null(start.time))
      times.vals <- times.vals[times.vals <= end.time]
    else
      if (is.null(end.time))
        times.vals <- times.vals[times.vals >= start.time]
      else
        times.vals <- times.vals[times.vals >= start.time & times.vals <= end.time]
      interval.resp <- getTimesSubset(response, data = data, which.times = times.vals, 
                                      times = times, include.times = TRUE)
      interval.resp <- cbind(getTimesSubset(individuals, data = data, which.times = times.vals, 
                                            times = times, 
                                            include.times = FALSE),
                             interval.resp)
  }
  
  #Calculate a value within an interval for each individual
  val.dat <- splitValueCalculate(response=response, weights=weights, individuals = individuals, 
                                 FUN = FUN, which.obs = which.obs, which.values = which.values, 
                                 data = interval.resp, na.rm=na.rm, sep=sep, ...)
  if (!is.null(suffix.interval))
  {
    names(val.dat)[match(paste(response, FUN, sep="."), names(val.dat))] <- 
      paste(response, FUN, suffix.interval, sep=".")
    if (which.obs)
      names(val.dat)[match(paste(response,FUN, "obs", sep="."), names(val.dat))] <- 
        paste(response, FUN, "obs", suffix.interval, sep=".")
    if (!is.null(which.values))
      names(val.dat)[match(paste(response, FUN, which.values, sep="."), names(val.dat))] <- 
        paste(response, FUN, which.values, suffix.interval, sep=".")
  }
  
  return(val.dat)
}

# Function to calculate water use indices (WUI) over an interval for a set of responses
"intervalWUI" <- function(responses, water.use = "Water.Use", individuals = "Snapshot.ID.Tag", 
                          times = "Days", start.time, end.time, suffix.interval = NULL, 
                          data, include.total.water = FALSE, na.rm = FALSE)
{ 
  .Deprecated(old = "intervalWUI", new = "byIndv4Intvl_WUI", package = "growthPheno",
              msg = paste0("'intervalValueCalculate' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'byIndv4Intvl_WUI' instead."))
  
  if (!all(sapply(list(start.time, end.time, suffix.interval), 
                  function(x) length(x) == 1)))
    stop("At least one of start.time, end.time and suffix.interval are not length one")
  
  #Check that response and individuals are in data
  if (!all(c(water.use, individuals, times) %in% names(data)))
    stop("One or more of water use, response and indivduals are not in data")
  
  #get the water use
  sum.dat <- splitValueCalculate(water.use, individuals = individuals, 
                                 FUN = "sum", na.rm = na.rm, 
                                 data = subset(data, 
                                               convertTimes2numeric(eval(parse(text=times))) >= 
                                                 max(start.time)+1 & 
                                               convertTimes2numeric(eval(parse(text=times))) <= 
                                                 max(end.time)))
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
  interval.resp <- cbind(getTimesSubset(individuals, data = data, times = times, 
                                        which.times = start.time, suffix = NULL),
                         getTimesSubset(responses, data = data, times = times, 
                                        which.times = start.time, include.times = TRUE, 
                                        suffix = "start"),
                         getTimesSubset(responses, data = data, times = times, 
                                        which.times = end.time, include.times = TRUE, 
                                        suffix = "end"))
  interval.resp <- merge(interval.resp, sum.dat, by = individuals)
  interval.resp <- interval.resp[ do.call(order, interval.resp), ]
  water.index <- lapply(responses, 
                        function(name, interval.resp, start.time, end.time, water.name, 
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
                        start.time = start.time, end.time = end.time, 
                        water.name = water.name, include.total.water = include.total.water)
  water.index <- as.data.frame(water.index)
  if (include.total.water)
    water.index <- cbind(interval.resp[c(individuals, water.name)], water.index)
  else 
    water.index <- cbind(interval.resp[individuals], water.index)
  return(water.index)
}

"plotTrait" <- function(tmp, response, response.smooth, x, xname, 
                        individuals, id.cols, times.factor, traits, 
                        methlabs, smethods, df, 
                        plots, devnplots, facet.x, facet.y, labeller = NULL, 
                        colour = "black", colour.column, colour.values = NULL, alpha, 
                        x.title, y.title = NULL, ggplotFuncs, ...)
{
  if (is.null(y.title))
    y.title <- response
  
  if (!("none" %in% plots))
  {
    #Plot comparisons
    if (any(c("methods+rawcompare", "methodscompare") %in% plots))
    {
      for (degfree in df)
      { 
        kresponses <- c(response, 
                        paste(response.smooth, methlabs[smethods], degfree, sep="."))
        tmp.sm <- tmp[c(id.cols, kresponses)]
        tmp.sm <- reshape(tmp.sm, direction = "long", 
                          varying = kresponses, v.names = response, 
                          idvar = id.cols, timevar = "Method")
        methods.df <- paste(methlabs, degfree, sep = "-")
        if ("methods+rawcompare" %in% plots)
        {
          if (length(smethods) == 1)
          {
            scale.labs <- c("Raw", methods.df)
          }
          else
          {
            if (length(smethods)  == 2)
            {
              scale.labs <- c(methods.df[1], "Raw", methods.df[2])
              tmp.sm <- within(tmp.sm, 
                               {
                                 Method <- fac.recode(factor(Method), c(2,1,3))
                                 Method <- factor(Method, labels = scale.labs)
                               })
            }
            else
              stop("A maximum of two smoothing.methods can be compared with raw data")
          }
          tmp.sm <- within(tmp.sm, 
                           {
                             Method <- factor(Method, labels = scale.labs)
                           })
          plt <- plotLongitudinal(data = tmp.sm, x=x, xname = xname, 
                                  response = response, individuals = individuals, 
                                  facet.x="Method", facet.y=facet.y, labeller = labeller, 
                                  colour = colour, colour.column = colour.column, 
                                  colour.values = colour.values, alpha = alpha, 
                                  x.title = x.title, y.title = y.title, 
                                  printPlot=FALSE, ggplotFuncs = ggplotFuncs, ...)
          print(plt)
        } else
        {
          tmp.sm <- within(tmp.sm, 
                           Method <- factor(Method, labels = c("Raw",  methods.df)))
          tmp.sm.sub <- tmp.sm[tmp.sm$Method != "Raw",]
          tmp.sm.sub <- within(tmp.sm.sub, 
                               Method <- factor(Method, labels = methods.df))
          plt <- plotLongitudinal(data = tmp.sm.sub, x=x, xname = xname, 
                                  response = response, individuals = individuals, 
                                  facet.x="Method", facet.y=facet.y, labeller = labeller, 
                                  colour = colour, colour.column = colour.column, 
                                  colour.values = colour.values, alpha = alpha, 
                                  x.title = x.title, y.title = y.title, 
                                  printPlot=FALSE, ggplotFuncs = ggplotFuncs, ...)
          print(plt)
        }
        
        #Plot deviation plots for both methods compare
        if (any(c("absolute.boxplots", "relative.boxplots") %in% devnplots))
        {
          y.titles <- c(paste("Absolute", response, "deviations for", degfree, "df", sep = " "),
                        paste("Relative", response, "deviations for", degfree, "df", sep = " "))
          names(y.titles ) <- c("absolute.boxplots", "relative.boxplots")
          y.titles <- y.titles[c("absolute.boxplots", "relative.boxplots") %in% devnplots]
          names(tmp.sm)[match(response, names(tmp.sm))] <- response.smooth
          tmp.sm <- merge(tmp.sm, tmp[c(id.cols, response)])
          tmp.sm <- subset(tmp.sm, Method != "Raw")
          plotDeviationsBoxes(data = tmp.sm, x.factor = times.factor, 
                              observed = response, smoothed = response.smooth, 
                              deviations.plots = devnplots, 
                              x.title = x.title, y.titles = y.titles, 
                              facet.x="Method", facet.y=facet.y, 
                              df = degfree)
        }
      } 
    } else
      if (any(c("df+rawcompare", "dfcompare") %in% plots))
      {
        for (smethod in smethods)
        {
          kresponses <- c(response, 
                          paste(response.smooth, methlabs[smethod], df, sep="."))
          tmp.sm <- tmp[c(id.cols, kresponses)]
          tmp.sm <- reshape(tmp.sm, direction = "long", 
                            varying = kresponses, v.names = response, 
                            idvar = id.cols, timevar = "DF")
          method.dfs <- paste(methlabs[smethod], df, sep = "-")
          if ("df+rawcompare" %in% plots)
          {
            if (length(df) != 2)
            {
              scale.labs <- c("Raw", method.dfs)
            }
            else
            {
              scale.labs <- c(method.dfs[1], "Raw", method.dfs[2])
              tmp.sm <- within(tmp.sm, 
                               {
                                 DF <- fac.recode(factor(DF), c(2,1,3))
                                 DF <- factor(DF, labels = scale.labs)
                               })
            }
            tmp.sm <- within(tmp.sm, 
                             {
                               DF <- factor(DF, labels = scale.labs)
                             })
            
            plt <- plotLongitudinal(data = tmp.sm, x=x, xname = xname, 
                                    response = response, individuals = individuals, 
                                    facet.x="DF", facet.y=facet.y, labeller = labeller, 
                                    colour = colour, colour.column = colour.column, 
                                    colour.values = colour.values, alpha = alpha, 
                                    x.title = x.title, y.title = y.title, 
                                    printPlot=FALSE, ggplotFuncs = ggplotFuncs, ...)
            print(plt)
          } else
          {
            tmp.sm <- within(tmp.sm, 
                             DF <- factor(DF, labels = c("Raw",  method.dfs)))
            tmp.sm.sub <- tmp.sm[tmp.sm$DF != "Raw",]
            tmp.sm.sub <- within(tmp.sm.sub, 
                                 DF <- factor(DF, labels = method.dfs))
            plt <- plotLongitudinal(data = tmp.sm.sub, x=x, xname = xname, 
                                    response = response, individuals = individuals, 
                                    facet.x="DF", facet.y=facet.y, labeller = labeller, 
                                    colour = colour, colour.column = colour.column, 
                                    colour.values = colour.values, alpha = alpha, 
                                    x.title = x.title, y.title = y.title, 
                                    printPlot=FALSE, ggplotFuncs = ggplotFuncs, ...)
            print(plt)
          }
          
          if (any(c("absolute.boxplots", "relative.boxplots") %in% devnplots))
          {
            #Plot deviation plots for both df compare
            y.titles <- c(paste("Absolute", response, "deviations for", smethod, "smoothing", sep = " "),
                          paste("Relative", response, "deviations for", smethod, "smoothing", sep = " "))
            names(y.titles ) <- c("absolute.boxplots", "relative.boxplots")
            y.titles <- y.titles[c("absolute.boxplots", "relative.boxplots") %in% devnplots]
            names(tmp.sm)[match(response, names(tmp.sm))] <- response.smooth
            tmp.sm <- merge(tmp.sm, tmp[c(id.cols, response)])
            tmp.sm <- tmp.sm[tmp.sm$DF != "Raw",]
            plotDeviationsBoxes(data = tmp.sm, x.factor = times.factor, 
                                observed = response, smoothed = response.smooth, 
                                deviations.plots = devnplots, 
                                x.title = x.title, y.titles = y.titles, 
                                facet.x="DF", facet.y=facet.y, 
                                df = NULL)
          }
          
        }
      } else
      {
        #Plot response
        if (("bothseparately" %in% plots) & ("response" %in% traits))
        { 
          pltu <- plotLongitudinal(data = tmp, x=x, xname = xname, 
                                   response = response, individuals = individuals, 
                                   facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                                   colour = colour, colour.column = colour.column, 
                                   colour.values = colour.values, alpha = alpha, 
                                   title="Unsmoothed response", x.title = x.title, y.title = y.title, 
                                   printPlot=FALSE, ggplotFuncs = ggplotFuncs, ...)
          print(pltu)
        }
        if (any(c("bothseparately", "smoothedonly") %in% plots))
        {
          #Plot smoothed response
          for (smethod in smethods)
          {
            for (degfree in df)
            { 
              r <- paste(response.smooth, methlabs[smethod], degfree, sep=".")
              plt <- plotLongitudinal(data = tmp, x=x, xname = xname, 
                                      response = r, individuals = individuals, 
                                      facet.x=facet.x, facet.y=facet.y, 
                                      labeller = labeller, 
                                      colour = colour, colour.column = colour.column, 
                                      colour.values = colour.values, alpha = alpha, 
                                      title="Smoothed response", x.title = x.title, y.title = r, 
                                      printPlot=FALSE, ggplotFuncs = ggplotFuncs, ...)
              print(plt)
              plotDeviationsBoxes(data = tmp, x.factor = times.factor, 
                                  observed = response, smoothed = r, 
                                  deviations.plots = devnplots, x.title = x.title, 
                                  facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                                  df = degfree)
            }
          }
        } else
        {
          if (any(c("absolute.boxplots", "relative.boxplots") %in% devnplots))
          {
            for (smethod in smethods)
            {
              for (degfree in df)
              { 
                #Plot deviation plots for both df compare
                y.titles <- c(paste("Absolute", response, "deviations for", 
                                    paste0(methlabs[smethod], degfree, sep="-"), sep = " "),
                              paste("Relative", response, "deviations for", 
                                    paste0(methlabs[smethod], degfree, sep="-"), sep = " "))
                names(y.titles ) <- c("absolute.boxplots", "relative.boxplots")
                y.titles <- y.titles[c("absolute.boxplots", "relative.boxplots") %in% devnplots]
                r <- paste(response.smooth, methlabs[smethod], degfree, sep=".")
                plotDeviationsBoxes(data = tmp, x.factor = times.factor, 
                                    observed = response, smoothed = r, 
                                    deviations.plots = devnplots, 
                                    x.title = x.title, y.titles = y.titles, 
                                    facet.x=facet.x, facet.y=facet.y, labeller = labeller, 
                                    df = degfree)
              }
            }
          }
          else
            stop(paste("which.plots option not allowed for:",plots))
        }
      }
  }
  invisible()
}

"probeSmoothing" <- function(data, response = "Area", response.smoothed = NULL, 
                             x = NULL, xname="xDays", 
                             times.factor = "Days", individuals="Snapshot.ID.Tag", 
                             na.x.action="exclude", na.y.action = "exclude", 
                             df, smoothing.methods = "direct", correctBoundaries = FALSE, 
                             get.rates = TRUE, rates.method="differences", 
                             facet.x = "Treatment.1", facet.y = "Smarthouse", 
                             labeller = NULL, x.title = NULL, 
                             colour = "black", colour.column=NULL, 
                             colour.values=NULL, alpha = 0.1, 
                             trait.types = c("response", "AGR", "RGR"), 
                             propn.types = c(0.1, 0.5, 0.75), propn.note = TRUE, 
                             which.plots = "smoothedonly",
                             deviations.plots = "none", alpha.med.devn = 0.5, 
                             ggplotFuncs = NULL, ggplotFuncsMedDevn = NULL, ...)
{ 
  .Deprecated(old = "probeSmoothing", new = "probeSmooths", package = "growthPheno",
              msg = paste0("'probeSmoothing' has been soft-deprecated; ",
                           "it is no longer being developed and will be removed in future versions.",
                           "\nUse 'probeSmooths' instead."))
  
  #check input arguments
  impArgs <- match.call()
  if ("na.rm" %in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  if ("smoothing.scales" %in% names(impArgs))
    stop("smoothing.scales has been deprecated; use smoothing.methods")
  if ("deviations.boxplots" %in% names(impArgs))
    stop("deviations.boxplots has been deprecated; use deviations.plots")
  options <- c("differences","derivative")
  opt <- options[check.arg.values(rates.method, options=options)]
  options <- c("none", "smoothedonly", "bothseparately", 
               "methodscompare", "methods+rawcompare", "dfcompare", "df+rawcompare")
  plots <- options[check.arg.values(which.plots, options=options)]
  if (any(c("bothseparately", "methodscompare", "dfcompare") %in% plots))
    plotunsmooth <- TRUE
  else
    plotunsmooth <- FALSE
  options <- c("none", "absolute.boxplots", "relative.boxplots", "compare.medians")
  devnplots <- options[unlist(lapply(deviations.plots, check.arg.values, 
                                     options=options))]
  if ("none" %in% devnplots & length(devnplots) > 1)
    devnplots <- "none"
  if (is.null(x.title))
    x.title <- times.factor
  
  options <- c("direct", "logarithmic")
  smethods <- options[unlist(lapply(smoothing.methods, check.arg.values, options=options))]
  methlabs <- c("Direct", "Log")
  names(methlabs) <- options
  methlabs <- methlabs[smethods]
  options <- c("response", "AGR", "RGR", "all")
  traits <- options[unlist(lapply(trait.types, check.arg.values, options=options))]
  if ("all" %in% traits)
    traits <- c("response", "AGR", "RGR")
  #If trait.types and get.rates don't match, go with option that is not the default
  grates <- NULL
  if (!any(c("AGR","RGR") %in% traits) && get.rates)
  {
    if (get.rates) 
    {
      get.rates <- FALSE
      warning("trait.types does not include AGR or RGR and so get.rates changed to FALSE")
    }
  } else
  {
    if (!get.rates)
    {
      if (length(traits) > 1 || traits != "response")
      {
        traits <- "response"
        warning("get.rates is FALSE and so trait.types changed to response")
      }
    }
    else
      grates <- c("AGR","RGR")[c("AGR","RGR") %in% traits]
  }
  #Form data.frame with just columns needed 
  id.cols <- c(individuals, times.factor, xname, colour.column)
  if (all(facet.x != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.x))
  if (all(facet.y != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.y))
  if (is.null(x))
    x <- xname
  v <- c(id.cols, response)
  #  else
  #    v <- c(v,x)
  if (!all(v %in% names(data)))
    stop(paste("Do not have the following required columns in data: ", 
               paste(v[!(v %in% names(data))],sep=", "), "\n", sep=""))
  tmp <- data[v]
  
  #Smooth response and form growth rates
  if (is.null(response.smoothed))
    response.smooth <- paste(response, "smooth", sep=".")
  else
    response.smooth <- response.smoothed
  responses.smooth <- response.smooth
  #no need for unsmoothed GRs if not plotted (need for dfcompare even though not plotted)
  #if ((plotunsmooth | !("none" %in% devnplots)) & get.rates) 
  #Always get rates if get.rates is TRUE so that they are in the returned data
  if (get.rates) 
  {
    tmp <- splitContGRdiff(tmp, response, individuals=individuals,
                           which.rates = grates, times.factor=times.factor, 
                           avail.times.diffs = FALSE)
  }
  for (degfree in df)
  { 
    for (smethod in smethods)
    {
      if (opt == "differences")
      { 
        tmp <- splitSplines(tmp, response, response.smoothed = response.smooth, 
                            x=xname, individuals = individuals, 
                            df = degfree, smoothing.method = smethod, 
                            correctBoundaries = correctBoundaries, 
                            na.x.action = na.x.action, na.y.action = na.y.action)
        if (get.rates)
        { 
          responses.smooth <- c(response.smooth, 
                                paste(response.smooth, grates, sep="."))
          tmp <- splitContGRdiff(tmp, response.smooth, individuals=individuals,
                                 which.rates = grates, times.factor=times.factor,
                                 avail.times.diffs = FALSE)
        } 
      } else #derivatives
      { 
        if (get.rates)
        { 
          responses.smooth <- c(response.smooth, 
                                paste(response.smooth, grates, sep="."))
          AGR <- NULL
          RGR <- NULL
          if ("AGR" %in% traits)
            AGR <- "AGR"
          if ("RGR" %in% traits)
            RGR <- "RGR"
          if (smethod == "direct")
            tmp <- splitSplines(tmp, response, response.smoothed = response.smooth, 
                                x=xname, individuals = individuals, deriv=1, 
                                suffices.deriv=AGR, extra.rate = "RGR", df = degfree, 
                                smoothing.method = smethod, 
                                na.x.action = na.x.action, na.y.action = na.y.action)
          else
            tmp <- splitSplines(tmp, response, response.smoothed = response.smooth, 
                                x=xname, individuals = individuals, deriv=1, 
                                suffices.deriv=RGR, extra.rate = "AGR", df = degfree, 
                                smoothing.method = smethod, 
                                na.x.action = na.x.action, na.y.action = na.y.action)
        }
        else
          tmp <- splitSplines(tmp, response, response.smoothed = response.smooth, 
                              x=xname, individuals = individuals, 
                              df = degfree, smoothing.method = smethod, 
                              correctBoundaries = correctBoundaries, 
                              na.x.action = na.x.action, na.y.action = na.y.action)
      }
      new.responses <- paste(responses.smooth, methlabs[smethod], degfree, sep=".")
      names(tmp)[match(responses.smooth, names(tmp))] <- new.responses
    }
  }
  
  #Plot some combination of unsmoothed and smoothed response, AGR and RGR
  if (!("none" %in% plots) | !("none" %in% devnplots))
  { 
    #Plot response
    if ("response" %in% traits)
      plotTrait(tmp = tmp, response = response, response.smooth = response.smooth, 
                x = x, xname = xname, individuals = individuals, 
                id.cols = id.cols, times.factor = times.factor, 
                traits = traits, methlabs = methlabs, smethods = smethods, df = df, 
                plots = plots, devnplots = devnplots, 
                facet.x = facet.x, facet.y = facet.y, labeller = labeller, 
                colour = colour, colour.column = colour.column, colour.values = colour.values, 
                alpha = alpha, x.title = x.title, y.title = response, 
                ggplotFuncs = ggplotFuncs, ...)
    
    #Plot growth rates
    for (grate in grates)
    {
      kresp <- paste(response, grate, sep = ".")
      kresp.sm <- paste(response.smooth, grate, sep = ".")
      plotTrait(tmp = tmp, response = kresp, response.smooth = kresp.sm, 
                x = x, xname = xname, individuals = individuals, 
                id.cols = id.cols, times.factor = times.factor, 
                traits = traits, methlabs = methlabs, smethods = smethods, df = df, 
                plots = plots, devnplots = devnplots, 
                facet.x = facet.x, facet.y = facet.y, labeller = labeller, 
                colour = colour, colour.column = colour.column, colour.values = colour.values, 
                alpha = alpha, x.title = x.title, y.title = kresp, 
                ggplotFuncs = ggplotFuncs, ...)
      
    }
  }
  
  if ("compare.medians" %in% devnplots)
  {
    plotMedianDeviations(tmp, response = response, response.smoothed = response.smooth, 
                         x = x, xname = xname, individuals = individuals, 
                         x.title = x.title, 
                         facet.x = facet.x, facet.y = facet.y, 
                         labeller = labeller, 
                         trait.types = traits, 
                         propn.types = propn.types, propn.note = propn.note, 
                         alpha.med.devn = alpha.med.devn,
                         smoothing.methods = smethods, df = df, 
                         ggplotFuncsMedDevn = ggplotFuncsMedDevn)  
    
  }
  invisible(tmp)
}
