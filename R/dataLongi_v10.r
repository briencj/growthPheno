globalVariables(c("Snapshot.ID.Tag", "Snapshot.Time.Stamp", "Time.after.Planting..d.", 
                  "Projected.Shoot.Area..pixels.", 
                  "Smarthouse", "DAP", "xDAP", "cPosn", "cMainPosn", "PSA", "Days",
                  "Genotype.ID", "Treatment.1", "Zone", "Lane", "ZLane",
                  "SHZone", "Mainunit", "ZMainunit", "cZone", 
                  "Hour", "Area.SV", "Area.SV1", "Area.SV2", "Area.TV", "Image.Biomass", "Centre.Mass", 
                  "Convex.Hull.SV", "Convex.Hull.TV", "Compactness.SV", "Max.Height", "Density", "Volume",
                  "Center.Of.Mass.Y.SV1", "Center.Of.Mass.Y.SV2", "Convex.Hull.Area.SV1",
                  "Convex.Hull.Area.SV1", "Convex.Hull.Area.TV", "Max.Dist.Above.Horizon.Line.SV1", 
                  "Max.Dist.Above.Horizon.Line.SV2", 
                  "Weight.Before", "Weight.After", "Water.Amount", "r", "Cumulative.Propn",
                  "Scheme", "DF"),
                "growthPheno", add = TRUE)

#Function to move camera prefix to a suffix without the characters up to the first _ (eg RBG_) 
#If no _ then moves whole prefix.
#Assumes that prefix ends at first full.stop
"pre2suffix" <- function(name, labsCamerasViews, keepCameraType)
{ 
  prefix <- (strsplit(name, ".", fixed=TRUE))[[1]][1]
  if (any(labsCamerasViews %in% prefix))
  { 
    if (grepl("_", prefix) && !keepCameraType)
    {
      suffix <- (strsplit(prefix, "_"))[[1]]
      if (length(suffix) == 2)
        suffix <- suffix[2]
      else
        suffix <- paste(suffix[2:length(suffix)], collapse = "_")
    }
    else
      suffix <- prefix
    fst <- nchar(prefix)+2
    name <- paste(substring(name, first=fst), suffix, sep=".")
  }
  return(name)
}

"importExcel" <- function(file, sheet = "raw data", sep = ",", cartId = "Snapshot.ID.Tag", 
                          imageTimes = "Snapshot.Time.Stamp", 
                          timeAfterStart = "Time.after.Planting..d.", 
                          cameraType = "RGB", keepCameraType = FALSE, 
                          labsCamerasViews = NULL, prefix2suffix = TRUE, 
                          startTime = NULL, timeFormat = "%Y-%m-%d %H:%M", 
                          plotImagetimes = TRUE, ...)
{ 
  #Check arguments
  impArgs <- match.call()
  if ("intervals" %in% names(impArgs))
    stop(paste("importExcel assumes that intervals are specified by timeAfterStart; \n", 
               "to have different intervals, call plotImagetimes separately"))
  if ("timeAfterPlanting"%in% names(impArgs))
    stop("timeAfterPlanting has been deprecated; use timeAfterStart")
  if ("planting.time"%in% names(impArgs))
    stop("planting.time has been deprecated; use startTime")
  
  #Input the raw imaging data
  if (grepl("csv", file))
  { 
    raw.dat <- read.csv(file, sep = sep, as.is=TRUE)
    raw.dat[imageTimes] <- as.POSIXct(raw.dat[[imageTimes]], format = timeFormat)
  }
  else if(grepl("xlsx", file))
  {
    raw.dat <- as.data.frame(read_excel(file, sheet=sheet))
    colnames(raw.dat) <- make.names(colnames(raw.dat))
    #raw.dat <- readWorksheetFromFile(file, sheet=sheet)
  }
  else
    stop("File name does not include csv or xlsx")
  ncinput <- ncol(raw.dat)
  if (!(cartId %in% names(raw.dat)))
    stop("cartId not in imported data")
  if (!(imageTimes%in% names(raw.dat)))
    stop("imageTimes not in imported data")
  
  
  #Rename image columns 
  #Change cameras and views as specified by labsCamerasViews
  vars <- names(raw.dat)
  if (!is.null(names(labsCamerasViews)))
  {
    old.names <- names(labsCamerasViews)
    for (old in old.names)
      vars <- gsub(old, labsCamerasViews[old], vars, fixed = TRUE)
    names(raw.dat) <- vars
  } else #find
  {
    labsCamerasViews <- vars[grep(cameraType, vars, fixed = TRUE)]
    if (length(labsCamerasViews) == 0)
      warning(paste("No imaging variables for a camera of type ", cameraType, " found", sep=""))
    labsCamerasViews <- strsplit(labsCamerasViews, ".", fixed=TRUE)
    labsCamerasViews <- unique(unlist(lapply(labsCamerasViews, 
                                             function(name) {return(name[[1]][1])})))
    names(labsCamerasViews) <- labsCamerasViews
  }
  
  #Move prefix for camera View to become a suffix without the RGB_
  if (prefix2suffix)
  { 
    newvars <- sapply(vars[1:length(vars)], pre2suffix, 
                      labsCamerasViews = labsCamerasViews, keepCameraType = keepCameraType)
    names(newvars) <- NULL
    names(raw.dat)[match(vars, names(raw.dat))] <- newvars
  } else
  {
    if (!keepCameraType)
    {
      vars <- names(raw.dat)
      vars <- gsub(paste(cameraType, "_",sep = ""), "", vars, fixed = TRUE) 
      names(raw.dat) <- vars
    }
  }
  
  #Change day calculation to take away a time origin and truncate to the nearest whole day
  #  if (!is.null(startTime))
  raw.dat <- calcTimes(raw.dat, imageTimes = imageTimes, timeFormat = timeFormat,
                       intervals = timeAfterStart , startTime = startTime,
                       intervalUnit = "days", timePositions = "Hour")

  #Plot the imaging times if required
  if (plotImagetimes)
    plotImagetimes(raw.dat, intervals=timeAfterStart, timePositions = "Hour", 
                   groupVariable = cartId, ...)
  
  #Check unique for Snapshot.ID.Tag, Time.after.Planting..d.
  combs <- as.vector(table(raw.dat[[cartId]], raw.dat[[timeAfterStart]]))
  if (any(combs != 1))
    warning(paste("There is not just one observation for",  
                  length(combs[combs != 1]), 
                  "combination(s) of",cartId,"and", timeAfterStart))
  
  #Sort data into cartId, Time.after.Planting..d. order and store
  # - may need to be reordered for analysis purposes
  raw.dat <- raw.dat[order(raw.dat[[cartId]], raw.dat[[timeAfterStart]]), ]
  return(raw.dat)
}


#Function to reduce imaging responses to those to be retained, forming image.dat
"prepImageData" <- function(data, cartId = "Snapshot.ID.Tag", 
                            imageTimes = "Snapshot.Time.Stamp", 
                            timeAfterStart = "Time.after.Planting..d.", 
                            PSAcolumn = "Projected.Shoot.Area..pixels.", 
                            idcolumns = c("Genotype.ID","Treatment.1"),
                            traits = list(all = c("Area", "Boundary.Points.To.Area.Ratio", 
                                                  "Caliper.Length", "Compactness", 
                                                  "Convex.Hull.Area"), 
                                          side = c("Center.Of.Mass.Y", 
                                                   "Max.Dist.Above.Horizon.Line")),
                            labsCamerasViews = list(all = c("SV1", "SV2", "TV"),
                                                    side = c("SV1", "SV2")), 
                            smarthouse.lev = NULL, 
                            calcWaterUse = TRUE)
{ 
  #Extract variables from data to form data frame of longitudinal data
  posndatevars <- c(cartId,timeAfterStart,
                    "Smarthouse","Lane","Position",imageTimes)
  imagevars <- NULL
  if (is.list(traits) && !is.null(traits))
  {
    if (is.null(labsCamerasViews))
      imagevars <- unlist(traits)
    else
    {
      if (!is.list(labsCamerasViews) || length(traits) != length(labsCamerasViews))
        stop(paste0("When traits is a list then labsCamerasViews must also be a list ",
                    "with the same number of components as traits"))
      if (length(labsCamerasViews) == 1)
        imagevars <- as.vector(outer(traits[[1]], labsCamerasViews[[1]], paste, sep = "."))
      else
        imagevars <- unlist(mapply(FUN =  function(traits, names) {
          if (is.null(names))
            t <- traits
          else
            t <- as.vector(t(outer(traits, names, paste, sep = ".")))
          invisible(t)
        }, 
        traits, labsCamerasViews))
      names(imagevars) <- NULL
    }
  } else
  {
    if (is.character(traits))
    {
      if (is.null(labsCamerasViews))
        imagevars <- unlist(traits)
      else
        imagevars <- as.vector(outer(traits, labsCamerasViews, paste, sep = "."))
    } else
      stop("traits is neither a list nor a character")
  }
  if (calcWaterUse)
    vars <- c(posndatevars, idcolumns, "Weight.Before","Weight.After","Water.Amount",
              PSAcolumn, imagevars)
  else
    vars <- c(posndatevars, idcolumns, PSAcolumn, imagevars)
  
  #Check that vars are in data
  if (!all(vars %in% names(data)))
    stop(paste("The following variables are not present in data:  ",
               paste(vars[!(vars %in% names(data))], collapse = ", "), sep = ""))
  
  image.dat <- data[, vars]
  
  #Add factors and variates needed in the analysis
  image.dat <- image.dat[do.call(order, image.dat), ]
  if (is.null(smarthouse.lev))
    smarthouse.lev <- unique(image.dat$Smarthouse)
  image.dat <- within(image.dat, 
                      { 
                        Smarthouse <- factor(Smarthouse, levels=smarthouse.lev)
                        xDAP <- as.numeric(image.dat[[timeAfterStart]])
                        DAP <- factor(xDAP, levels = sort(unique(xDAP)))
                        cPosn <- Position - mean(unique(Position))
                        Position <- factor(Position, levels=sort(unique(Position)))
                        PSA <- image.dat[[PSAcolumn]]/1000
                      })
  
  facs <- c("Lane", idcolumns, "DAP")
  image.dat[facs] <- as.data.frame(lapply(image.dat[facs], FUN = factor))
  
  
  #Now derive a Reps factor 
  #+
  if (all(idcolumns %in% vars))
  {
    image.dat <- within(image.dat, 
                        { 
                          Reps <- 1
                          trts <- dae::fac.combine(as.list(image.dat[idcolumns]))
                        })
    for (t in levels(image.dat$trts))
    { 
      which.indiv <- with(image.dat, 
                          sort(unique(image.dat[trts==t, cartId])))
      for (k in 1:length(which.indiv))
        image.dat[image.dat$trts == t & 
                    image.dat$Snapshot.ID.Tag == which.indiv[k], "Reps"] <- k
    }
    image.dat$Reps <- factor(image.dat$Reps)
  } else 
    image.dat$Reps <- NA
  
  #Form responses that can be calculated by row-wise  operations: 
  image.dat <- calcTimes(image.dat, imageTimes = imageTimes,
                         timePositions = "Hour")
  kpx.vars <- imagevars[c(grep("Area.", imagevars, fixed = TRUE), 
                          grep("Convex.Hull.Circumference", imagevars, fixed = TRUE))]
  kpx.vars <- kpx.vars[!grepl("Ratio", kpx.vars, fixed = TRUE)]
  image.dat[kpx.vars] <- image.dat[kpx.vars]/1000
  
  #'## Calculate Water Use
  #+
  if (calcWaterUse)
    image.dat <- within(image.dat, 
                        { 
                          WU <- unlist(by(Weight.After, list(Snapshot.ID.Tag), 
                                          FUN=calcLagged)) - Weight.Before
                        })
  
  
  out.posndatevars <- c("Smarthouse","Lane","Position","DAP", "xDAP",
                        cartId, imageTimes, "Hour", "Reps")
  imagevars <- c("PSA", imagevars)
  if (calcWaterUse)
    imagevars <- c("Weight.Before","Weight.After","Water.Amount", "WU", imagevars)
  out.vars <- c(out.posndatevars, idcolumns, imagevars)
  
  #Re-order rows and response columns
  image.dat <- image.dat[order(image.dat[[cartId]], image.dat$DAP), ]
  image.dat <- image.dat[out.vars]
  names(image.dat) <- gsub("Area", "PSA", names(image.dat), fixed = TRUE)
  return(image.dat)
}

#Function to add design factors for blocked, split plot design
"designFactors" <- function(data, insertName = NULL, designfactorMethod = "LanePosition", 
                            nzones = 6, nlanesperzone = 4, nmainunitsperlane = 11, nsubunitspermain = 2)
{ 
  options <- c("LanePosition","StandardOrder")
  desfactor <- options[check.arg.values(designfactorMethod, options=options)]
  
  #Extract variables from data
  vars <- names(data)
  required <- c("Smarthouse", "Snapshot.ID.Tag", "xDAP")
  if (desfactor == "LanePosition")
    required <- c("Lane", "Position", required)
  checkNamesInData(required, data)
 
  n <- nrow(data)
  smarthouse.lev <- levels(data$Smarthouse)
  nshouse <- length(smarthouse.lev)
  ncarts = length(unique(data$Snapshot.ID.Tag))
  nexpcarts <- nshouse * nzones * nlanesperzone * nmainunitsperlane * nsubunitspermain
  
  if (desfactor == "StandardOrder" )
  { if (ncarts != nexpcarts)
    stop(paste("The number of unique Snapshot.ID.Tag values must be equal to the product of the numbers of\n",
               "      smarthouses, Zone, lanes per zone, mainunits per zone and subplots per mainunit"))
  } else
  { if (ncarts != nexpcarts)
    warning(paste("The number of unique Snapshot.ID.Tag values is not equal to the product of the numbers of\n",
                  "      smarthouses, Zone, lanes per zone, mainunits per zone and subplots per mainunit"))
  }
  if (n %% ncarts != 0)
    warning("There is not the same number imagings for each cart")
  
  #Add factors and variates needed in the analysis
  data <- data[do.call(order, data), ]
  
  #Generate design factors
  if (desfactor == "LanePosition")
  { if (!is.factor(data$Lane))
  {
    levs <- unique(data$Lane)
    levs <- levs[order(levs)]
    data <- cbind(data, 
                  with(data, fac.divide(factor(Lane, levels = levs), 
                                        list(Zone=nzones, ZLane = nlanesperzone))))
  } else
    data <- cbind(data, 
                  with(data, fac.divide(Lane, 
                                        list(Zone=nzones, ZLane = nlanesperzone))))
  if (!is.factor(data$Position))
  {
    levs <- unique(data$Position)
    levs <- levs[order(levs)]
    data <- cbind(data, 
                  with(data, 
                       fac.divide(factor(Position, levels = levs), 
                                  list(Mainunit=nmainunitsperlane, 
                                       Subunit = nsubunitspermain))))
  } else
    data <- cbind(data, 
                  with(data, 
                       fac.divide(Position, 
                                  list(MainUnit=nmainunitsperlane, 
                                       Subunit = nsubunitspermain))))
  } else
    if (desfactor == "StandardOrder")
    { id <- unique(data$Snapshot.ID.Tag)
    data <- merge(data, 
                  data.frame(fac.gen(list(Smarthouse=smarthouse.lev, 
                                          Zone=nzones, ZLane = nlanesperzone, 
                                          Mainunit=nmainunitsperlane, 
                                          Subunit = nsubunitspermain))[,-1],
                             Snapshot.ID.Tag = id), 
                  all.x=TRUE, sort=FALSE)
    } 
  data <- within(data, 
                 {
                   cPosn <- dae::as.numfac(Position)
                   cPosn <- cPosn - mean(unique(cPosn))
                 })
  if (nshouse == 1)
  {
    xMain <- with(data, aggregate(cPosn, by=list(Zone, ZLane, Mainunit), mean))
    names(xMain) <- c("Zone", "ZLane", "Mainunit", "cMainPosn") 
    data <- merge(data, xMain, all.x =TRUE, by = c("Zone", "ZLane", "Mainunit"), sort=FALSE)
    
  } else
  {
    xMain <- with(data, aggregate(cPosn, by=list(Smarthouse, Zone, ZLane, Mainunit), mean))
    names(xMain) <- c("Smarthouse", "Zone", "ZLane", "Mainunit", "cMainPosn") 
    data <- merge(data, xMain, all.x =TRUE, sort=FALSE)
  }
  data <- with(data, data[order(Snapshot.ID.Tag, xDAP), ])
  data <- within(data, { SHZone <- fac.combine(list(Smarthouse,Zone))
  ZMainunit <- fac.combine(list(ZLane,Mainunit))
  cZone <- as.numeric(Zone)
  cZone <- cZone - mean(unique(cZone))
  
  })
  
  
  facs <- c("Zone","cZone","SHZone","ZLane","ZMainunit",
            "Subunit", "cMainPosn","cPosn")
  out.vars <- c(vars,facs)
  if (!is.null(insertName))
  { k <- match(insertName, vars)
  if (!is.na(k))
    out.vars <- c(vars[1:k],facs,vars[(k+1):length(vars)])
  }
  
  #Re-order rows and response columns
  data <- with(data, data[order(Snapshot.ID.Tag, xDAP), ])
  data <- data[out.vars]
  return(data)
}


#Function that produces a longitudinal plot
"plotProfiles" <- function(data, response = "PSA", 
                           individuals="Snapshot.ID.Tag", times = "DAP", 
                           x = NULL, title = NULL, 
                           x.title = "DAP", y.title = "PSA (kpixels)", 
                           facet.x = ".", facet.y =   ".", 
                           labeller = NULL, scales = "fixed", 
                           breaks.spacing.x = -2, angle.x = 0, 
                           colour = "black", 
                           colour.column=NULL, colour.values=NULL, 
                           alpha = 0.1, addMediansWhiskers = FALSE, 
                           ggplotFuncs = NULL, 
                           printPlot = TRUE)
{ 
  strip.text.size <- 10
  
  checkNamesInData(c(response,individuals,times),data)

  data <- data[!is.na(data[response]),]
  data[times] <- convertTimes2numeric(data[[times]])
  if (is.null(x))
    x <- times

  longi.plot <- ggplot(data=data, aes(x = .data[[!!x]], y = .data[[!!response]])) +
    theme_bw() +
    setScaleTime(data[[x]], breaks.spacing.x = breaks.spacing.x) +
    theme(panel.grid.major = element_line(colour = "grey60", linewidth = 0.5), 
          panel.grid.minor = element_line(colour = "grey80", linewidth = 0.5),
          axis.text.x = element_text(size = 7.5, angle = angle.x)) +
    xlab(x.title) + ylab(y.title) + ggtitle(title)
  
  #Do facet if have any
  if (all(facet.x != ".") | all(facet.y != "."))
  {
    facet.form <- facet.char2formula(facet.x, facet.y)
    if (is.null(labeller))
      longi.plot <- longi.plot + facet_grid(facet.form, scales = scales)
    else
      longi.plot <- longi.plot + facet_grid(facet.form, labeller = labeller, scales = scales)
    longi.plot <- longi.plot + theme(strip.text = element_text(size=strip.text.size, face="bold"),
                                     axis.title = element_text(face="bold"), legend.position="none")
  }
  if (is.null(colour.column))
    longi.plot <- longi.plot + geom_line(aes(group=.data[[!!individuals]]),  
                                         colour=colour, alpha=alpha)
  else
    longi.plot <- longi.plot + geom_line(aes(group=.data[[!!individuals]], 
                                             colour=.data[[!!colour.column]]), 
                                         alpha=alpha)
  if (!(is.null(colour.values)))
    longi.plot <- longi.plot + scale_colour_manual(values = colour.values)
  
  
  #Calculate the medians and outer whisker-values over time and facetting factors
  if (addMediansWhiskers)
  {
    if (x != times)
      warning("x is ", x, " and times is ", times, "\nIs the value of times the name of the column from which x is derived?")
    
    #Create a factor Times that has the plotted values of x for its labels
    times.factor <- "Times"
    data[times.factor] <- data[times]
    data[times.factor] <- factor(unlist(data[times.factor]), 
                                 labels = unique(data[times.factor])[
                                   order(unique(data[[times.factor]])),])
    
    #Get facet cols if have any
    facet.cols <- NULL
    if (all(facet.x != ".") | all(facet.y != "."))
    {
      facet.cols <- c(facet.x, facet.y)
      facet.cols <- facet.cols[facet.cols != "."]
    }
    
    stats <- c("median", "lower.whisker", "upper.whisker")
    
    dat.split <- split(data, f = as.list(data[c(facet.cols, times.factor)]),
                       lex.order = TRUE)
    summ <- lapply(dat.split, 
                   function(x, response)
                   { 
                     stats <- boxplot.stats(as.vector(x[[response]]))$stats[c(3, 1, 5)]
                     return(stats)
                   },
                   response = response)
    summ <- as.data.frame(do.call(rbind, summ))
    summ <- cbind(fac.gen(lapply(data[c(facet.cols, times.factor)], levels)), 
                  summ)
    names(summ)[(ncol(summ)-2):ncol(summ)] <- stats
    summ[times] <- dae::as.numfac(unlist(summ[times.factor]))
    
    summ <- reshape(summ, direction = "long", 
                    varying = stats, 
                    v.names = response, 
                    idvar = c(facet.cols, times.factor), timevar = "Statistic")
    summ$Statistic <- factor(summ$Statistic, labels = stats)
    longi.plot <- longi.plot + 
      geom_line(data = summ[summ$Statistic == "median", ], 
                aes(x=.data[[x]], y=.data[[response]], alpha = 0.75),
                show.legend = FALSE, linetype="solid", colour = "black") +
      geom_line(data = summ[summ$Statistic != "median", ], 
                aes(x=.data[[x]], y=.data[[response]], group=.data[["Statistic"]], 
                    alpha = 0.75),
                show.legend = FALSE, linetype="dashed", colour = "black")
  }
  
  if (!is.null(ggplotFuncs))
    for (f in ggplotFuncs)
      longi.plot <- longi.plot + f
  
  if (printPlot)
    print(longi.plot)
  invisible(longi.plot)
}


#Function that calculates intervals and imageTimes from imageTimes
"calcTimes" <- function(data, imageTimes = NULL, 
                        timeFormat = "%Y-%m-%d %H:%M",
                        intervals = "Time.after.Planting..d.", startTime = NULL, 
                        intervalUnit = "days", timePositions = NULL)
{
  if (!is.null(imageTimes))
  {
    if (!(imageTimes %in% names(data)))
      stop("A column for imageTimes is not present in data")
    if (any(class(data[[imageTimes]])[1] %in% c("character", "factor")))
      data[imageTimes] <- as.POSIXct(data[[imageTimes]], format = timeFormat)
    units <- c("secs", "mins", "hours", "days")
    unit <- units[check.arg.values(intervalUnit, options=units)]
    if (unit == "secs")
    {
      d <- getOption("digits.secs")
      if (d == 0)
        warning(paste("Fractions of sections will not be stored or extracted unless: \n",
                      "(i) option(digits.secs) has been set to the number of decimal places required \n",
                      "and (ii) %OS is used for seconds in timeFormat",
                      sep=""))
    }
    if (!is.null(startTime))
    { 
      startTime <- as.POSIXct(startTime, format = timeFormat, tz = "UTC")
      data[[intervals]] <- difftime(data[[imageTimes]], startTime, units=intervalUnit)
      data[[intervals]] <- as.numeric(trunc(data[[intervals]], units=intervalUnit))
    }
    if (!is.null(timePositions))
    {
      data[[timePositions]] <- trunc(data[[imageTimes]], units=unit)
      if (unit == "secs")
      {
        data[[timePositions]] <- as.numeric(format(data[[imageTimes]], "%OS"))
        data[[timePositions]] <- data[[timePositions]] - floor(data[[timePositions]])
      }
      else
        data[[timePositions]] <- as.numeric(difftime(data[[imageTimes]], 
                                                     data[[timePositions]], 
                                                     units=units[(match(unit, units) - 1)]))
    }
  }
  return(data)
}

#Function that produces a plot of the imaging times
"plotImagetimes" <- function(data, intervals = "Time.after.Planting..d.", 
                             timePositions = "Hour", 
                             groupVariable = "Snapshot.ID.Tag", colourVariable = "Lane", 
                             ggplotFuncs = NULL, printPlot = TRUE)
{ 
  #Check whether have enough information to do the calculations
  if (!all(c(intervals, timePositions, groupVariable, colourVariable) %in% names(data)))
  {
    miss.names <- c(intervals, timePositions, groupVariable, 
                    colourVariable)[!(c(intervals, timePositions, groupVariable, colourVariable) 
                                      %in% names(data))]
    stop(paste(paste(miss.names, collapse = ", ")," is/are not present in the data", sep=""))
  }
   if (!(is.numeric(data[[intervals]])))
    data[intervals] <- dae::as.numfac(data[[intervals]])
  if (!(is.numeric(data[[colourVariable]])))
    data[colourVariable] <- dae::as.numfac(data[[colourVariable]])
  
  #Do plot
  start <- min(data[intervals], na.rm=TRUE)
  end <- max(data[intervals], na.rm=TRUE)
  time.plot <- ggplot(data, aes_string(x=intervals, y=timePositions)) +
    geom_line(aes_string(group=groupVariable, colour=colourVariable), alpha=0.05) + 
    scale_colour_gradient(low="grey60", high="grey20") + 
    geom_point(aes_string(group=groupVariable), size=0.5) +
    facet_grid(Smarthouse ~ .) + theme_bw() +
    scale_x_continuous(breaks=seq(start, end, by=2)) +
    ylab("Hour of day")
  
  if (!is.null(ggplotFuncs))
    for (f in ggplotFuncs)
      time.plot <- time.plot + f
  
  if (printPlot)
    print(time.plot)
  invisible(time.plot)
}

#New function
"plotAnom" <- function(data, response = "sPSA", 
                       individuals="Snapshot.ID.Tag", 
                       times = "DAP", x = NULL, 
                       breaks.spacing.x = -2, angle.x = 0, 
                       vertical.line=NULL, 
                       groupsFactor=NULL, lower=NULL, upper=NULL, 
                       start.time=NULL, end.time=NULL,  
                       suffix.interval=NULL, 
                       columns.retained=c("Snapshot.ID.Tag", "Smarthouse", "Lane", "Position", 
                                          "Treatment.1", "Genotype.ID"),
                       whichPrint=c("anomalous","innerPlot","outerPlot"), na.rm=TRUE, ...)
{ 
  inargs <- list(...)
  if ("breaks" %in% names(inargs))
    stop("The argument breaks has been replaced by breaks.spacing.x")

  if (!all(individuals %in% columns.retained))
    stop("The individuals column(s) is (are) not in the columns.retained")
  if (is.null(lower) & is.null(upper))
    stop("Must set at least one of lower and upper")
  options <- c("anomalous","innerPlot","outerPlot")
  opt <- options[unlist(lapply(whichPrint, check.arg.values, options=options))]
  
  #Set x and make times numeric
  if (is.null(x))
    x <- times
  data[times] <- convertTimes2numeric(data[[times]])
  
  #Determine anomalous individuals
  if (is.null(groupsFactor))
  { 
    anomalous.individuals <- byIndv4Intvl_ValueCalc(data = data, response = response, 
                                                    individuals = individuals, times = times, 
                                                    FUN = "anom",  
                                                    lower = lower, upper = upper, 
                                                    start.time = start.time, end.time = end.time, 
                                                    suffix.interval = suffix.interval, 
                                                    na.rm = na.rm)
    data <- merge(data, anomalous.individuals[,1:2], by=individuals, sort=FALSE)
  }
  else
  { 
    tmp <- split(data, data[[groupsFactor]])
    ngrps <- length(tmp)
    nstart <- length(start.time)
    nend <- length(end.time)
    nlow <- length(lower)
    nup <- length(upper)
    if (nstart == 0)
    { 
      kstart <- NULL
      if (nend != 1 & nend != ngrps)
        stop("Number of end.time values must be equal 1 or the number of levels in groupsFactor")
      kend <- end.time[1]
    } else
      if (nend ==0)
      { 
        kend <- NULL
        if (nstart != 1 & nstart != ngrps)
          stop("Number of start.time values must be equal 1 or the number of levels in groupsFactor")
        kstart <- start.time[1]
      } else
      {     
        if (nstart != nend | (nstart != ngrps & nstart != 1))
          stop("Number of start.time and end.time values must be equal and equal to 1 \n",
               "or the number of levels in groupsFactor")
        kstart <- start.time[1]
        kend <- end.time[1]
      }
    if (!(nlow == 0 | nlow  == ngrps |  nlow == 1))
      stop("Number of lower values must equal to 1 or the number of levels in groupsFactor")
    if (!(nup == 0 | nup  == ngrps | nup == 1))
      stop("Number of upper values must equal to 1 or the number of levels in groupsFactor")
    klow <- lower[1]
    kup <- upper[1]
    for (k in 1:ngrps)
    { 
      if (nstart > 1)
        kstart <- start.time[k]
      if (nend > 1)
        kend <- end.time[k]
      if (nlow > 1)
        klow <- lower[k]
      if (nup > 1)
        kup <- upper[k]
      anomalous.individuals <- byIndv4Intvl_ValueCalc(data = tmp[[k]], response = response, 
                                                      individuals = individuals, times = times, 
                                                      FUN = "anom", 
                                                      lower = klow, upper = kup, 
                                                      start.time = kstart, end.time = kend, 
                                                      suffix.interval = suffix.interval, 
                                                      na.rm = na.rm)
      tmp[[k]] <- merge(tmp[[k]], anomalous.individuals[,1:2], by=individuals, sort=FALSE)
    }
    data <- do.call(rbind, tmp)
  }
  response.anom <- names(anomalous.individuals)[2]
  
  #Plot without anomalous individuals
  if (sum(!na.omit(data[[response.anom]]) > 0))
  { 
    innerPlot <- plotProfiles(data = subset(data, !na.omit(data[[response.anom]])), 
                              response = response, times = times, x=x, 
                              breaks.spacing.x = breaks.spacing.x, 
                              printPlot=FALSE, ...)
    if (!is.null(vertical.line))
      innerPlot <- innerPlot + geom_vline(xintercept=vertical.line, linetype="longdash", linewidth=1)
    if ("innerPlot" %in% opt)
      print(innerPlot)
  } else
    innerPlot <- NULL
  
  #Print out anomalous individuals
  if ("anomalous" %in% opt)
  { 
    anom.dat <- data[c(columns.retained, response.anom)] 
    anom.dat <- split(anom.dat, anom.dat[[individuals]])
    anom.dat <- lapply(anom.dat, 
                       function(dat)
                         dat <- dat[1,])
    anom.dat <- do.call(rbind, anom.dat)
    anom.dat <- anom.dat[anom.dat[[response.anom]],]
    anom.dat <- anom.dat[order(anom.dat[[individuals]]), columns.retained]
    print(anom.dat)
  }  
  
  #Plot anomalous individuals, adding Snapshot.ID.Tag
  if (sum(data[[response.anom]] > 0, na.rm = TRUE))
  { 
    outerPlot <- plotProfiles(data = subset(data, data[[response.anom]]), 
                              times = times, x=x, response = response, alpha=0.5, colour="purple", 
                              breaks.spacing.x = breaks.spacing.x, 
                              printPlot=FALSE, ...)
    
    if (!is.null(vertical.line))
      outerPlot <- outerPlot + geom_vline(xintercept=vertical.line, linetype="longdash", linewidth=1)
    
    if ("outerPlot" %in% opt)
      print(outerPlot)
  } else
    outerPlot <- NULL
  
  invisible(list(data = data, innerPlot = innerPlot, outerPlot = outerPlot))
}

