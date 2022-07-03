#Function to check that arguments that come in via the ellipsis are all arguments to the allowed functions
#Usually funcs should be the current function and any to which ... is allowed to pass arguments
checkEllipsisArgs <- function(funcs, inargs)
{
  inargs <- names(inargs)
  if (length(inargs))
  {
    args <- unique(unlist(lapply(funcs, formalArgs)))
    args <- args[-match("...", args)]
    foreignArgs <- inargs[!(inargs %in% args)]
    if (length(foreignArgs))
      stop("the argument(s) ", paste0(foreignArgs, collapse = ", "), " are not legal arguments for ", 
           paste0(paste0("'",funcs, collapse = "', "), "'"))
  }
  invisible()
}

#Function to get the number of combinations of the factors in cols that are in data
getNschemes <- function(data, cols)
{
  if (!all(is.factor(data[cols])))
    data[cols] <- lapply(data[cols], as.factor)
  nsch <- length(levels(dae::fac.combine(as.list(data[cols]), 
                                         combine.levels = TRUE, sep = "-")))
}

#Function to have commas between non-smoothing-parameter factors and 
#hyphens between smoothing-parameter factors 
fac.mixcombine <- function(dat, factors, smooth.cols)
{
  tune.facs <- factors[factors %in% smooth.cols]
  if ("Method" %in% tune.facs)
    dat["Method"] <- dae::fac.recast(dat[["Method"]], 
                                     newlevels = substring(levels(dat[["Method"]]), 
                                                           first = 1, last = 3))
  if (length(tune.facs) > 0)
  {
    factors <- setdiff(factors, tune.facs)
    fac <- with(dat, dae::fac.combine(as.list(dat[tune.facs]),
                                      combine.levels = TRUE, sep = "-"))
    dat$tune.fac <- fac
    if (length(factors) > 0)
      fac <- with(dat, dae::fac.combine(as.list(dat[c(factors,"tune.fac")]),
                                        combine.levels = TRUE, sep = ", "))
  }
  else
    fac <- with(dat, dae::fac.combine(as.list(dat[factors]),
                                      combine.levels = TRUE, sep = ", "))
  return(fac)
}

#Functions for setting up and checking a smooths.frame

"is.smooths.frame" <- function(object)
{
  inherits(object, "smooths.frame") && inherits(object, "data.frame")
}


"validSmoothsFrame" <- function(object)
{
  smooth.cols <- c("Type","TunePar","TuneVal","Tuning","Method")
  issmoothsframe <- TRUE 
  #Check have only legal attributes
  if (!all(names(attributes(object)) %in% c("names", "row.names", "class", "n", "t", "nschemes", 
                                            "individuals", "times")))
  {
    issmoothsframe[1] <- FALSE
    issmoothsframe <- c(issmoothsframe, 
                        "\n  An unexpected attribute is present in a smooths.frame")
  }
  #Check that have non-null values for smooths.frames attributes
  smooth.attr <- c("n", "t", "nschemes", "individuals", "times")
  which.null <- unlist(lapply(smooth.attr, function(x, object) is.null(attr(object, which = x)), 
                              object = object))
  if (any(which.null))
  {
    issmoothsframe[1] <- FALSE
    issmoothsframe <- c(issmoothsframe, 
                        paste0("\n  The followng attributes of a smooths.frame are NULL: ", 
                               paste(smooth.attr[which.null],collapse=", ")))
  }
  #Check that is a data.frame
  if (!is.data.frame(object))
  {
    issmoothsframe[1] <- FALSE
    issmoothsframe <- c(issmoothsframe, 
                        "\n  smooths.frame is not a data.frame")
  }
  #Check have all the smoothing parameter columns 
  if (!all(smooth.cols %in% names(object)))
  {
    issmoothsframe[1] <- FALSE
    issmoothsframe <- c(issmoothsframe, 
                        paste0("\n  Do not have the following required smoothing-parameters columns in a smooths.frame: ", 
                               paste(smooth.cols[!(smooth.cols %in% names(object))],collapse=", ")))
  } else
  {
    #Check that the number of smoothing parameter sets in object is equal to the nschemes attribute
    nschemes <- getNschemes(object, smooth.cols)
    if (attr(object, which = "nschemes") != nschemes)
    {
      issmoothsframe[1] <- FALSE
      issmoothsframe <- c(issmoothsframe, 
                          paste0("\n The number of different combinations of the smoothing-parameter values ", 
                                 "in the smooths.frame does not match the number in its 'nschemes' attribute"))
    }
  }
  
  #Check the times.cols
  id.cols <- c(attr(object, which = "times"), attr(object, which = "individuals"))
  if (!(all(id.cols %in% names(object))))
  { 
    issmoothsframe[1] <- FALSE
    issmoothsframe <- c(issmoothsframe, 
                        paste0("\n  Do not have columns  for ", 
                               paste(id.cols[!(id.cols %in% names(object))],collapse=" or "), 
                               " in the smooths.frame"))
  }
  
  if (length(issmoothsframe) > 1)
    issmoothsframe[1] <- "Error in validSmoothsFrame : "
  return(issmoothsframe)
}

"as.smooths.frame" <- function(data, individuals = NULL, times = NULL)
{
  smooth.cols <- c("Type","TunePar","TuneVal","Tuning","Method")
  if (!all(smooth.cols %in% names(data)))
    stop(paste0("Cannot assign smooths.frame class to supplied data.frame",  
                " because it does not contain the following smoothing-parameters columns: ", 
                paste0(smooth.cols[!(smooth.cols %in% names(data))],collapse=", ")))
  #Set attributes of data
  class(data) <- c("smooths.frame", class(data))
  attr(data, which = "individuals") <- individuals
  attr(data, which = "n") <- length(unique(data[[individuals]]))
  attr(data, which = "times") <- times
  attr(data, which = "t") <- length(unique(data[[times]]))
  attr(data, which = "nschemes") <- getNschemes(data, smooth.cols)
  
  return(data) 
}  

methods::setOldClass("smooths.frame")


checkLayoutArgs <- function(data, plts.by, plts.group, facet.x, facet.y)
{  
  if (any(facet.x == ".")) facet.x <- NULL
  if (any(facet.y == ".")) facet.y <- NULL
  
  subfacs <- c(plts.by, plts.group, facet.x, facet.y)
  
  if (!is.allnull(subfacs) && !is.null(data) && !all(unique(subfacs) %in% names(data)))
    stop("The factor(s) ", paste0(unique(subfacs)[!(unique(subfacs) %in% names(data))], collapse = ", "), 
         " are not included in the smooths.frame")
  
  norepeats <- TRUE
  if (length(unique(subfacs)) < length(subfacs))
  {
    repfacs <- names(table(subfacs)[table(subfacs) > 1])
    stop("The factor(s) ", paste0(repfacs, collapse = ", "), 
         " occur(s) in more than one of the facet/plots arguments")
  }
  return(norepeats)
}

checkPlotsArgs <- function(data, plts.by, plts.group = NULL, facet.x, facet.y)
{
  smooth.cols <- c("Type","TunePar","TuneVal","Tuning","Method")
  
  #Check that smoothing-parameter factors in the plots subsetting arguments account for 
  # all of the smoothing parameter sets in data
  subfacs <- unique(c(plts.by, plts.group, facet.x, facet.y))
  if (is.allnull(subfacs))
    stop("There are no factors assigned to the plots and facet arguments ", 
         "- enough smoothing-parameter factors need to be assigned so that they uniquely index the combinations ", 
         "of the smoothing-parameter values in the smooths.frame")
  else
  {
    spar.facs <- subfacs[subfacs %in% smooth.cols]
    if (!any(smooth.cols %in% subfacs))
      stop("There are no smoothing-parameter factors assigned to the plots and facet arguments ", 
           "- enough of them need to be assigned so that they uniquely index the combinations ", 
           "of the smoothing-parameter values in the smooths.frame")
  }
  
  ncombos <- getNschemes(data, spar.facs)
  if (ncombos != attr(data, which = "nschemes"))
    stop(paste0("\n The number of different combinations of the smoothing-parameter values in ", 
                "the smooths.frame and of the levels combination of the following factors nominated ",
                "in the facet/plots arguments are not equal: ", paste0(spar.facs, collapse = ", ")))
  return(data)
}  

#Function to produce a single plot of deviations boxplots 
plotDeviationsBoxes <- function(data, observed, smoothed, x.factor, 
                                x.title = NULL, y.titles = NULL,
                                facet.x = ".", facet.y = ".", labeller = NULL, 
                                df, deviations.plots = "absolute", 
                                ggplotFuncs = NULL, printPlot = TRUE, ...)  
{
  options <- c("none", "profiles", "absolute.boxplots", "relative.boxplots", 
               "medians.deviations", "compare.medians") #needed for probeSmoothing only
  devnplots <- options[unlist(lapply(deviations.plots, check.arg.values, 
                                     options=options))]
  if (length(df) > 1)
    stop("Deviations boxplots can only be plotted for a single df at a time")
  if ("none" %in% devnplots)
    devnplots <- "none"
  else
  {
    devnplots <- devnplots[grepl("boxplots", devnplots, fixed = TRUE)]
    if (length(devnplots) == 0)
      devnplots <- "none"
  }
  
  plts <- list()
  #only do deviations boxplots if  requested  
  if (all(devnplots != "none"))
  {
    dat <- data
    if (!all(c(observed, smoothed) %in% names(dat)))
      stop(paste("One or more of", observed, "and", smoothed, "is missing from ", 
                 deparse(substitute(data))))
    if (is.null(x.title))
      x.title <- x.factor
    strip.text.size <- 10
    
    if (is.null(y.titles) && !("none" %in% devnplots))
    {
      y.titles <- c(paste("Absolute", observed, "deviations", sep = " "),
                    paste("Relative", observed, "deviations", sep = " "))
      names(y.titles ) <- c("absolute.boxplots", "relative.boxplots")
      y.titles <- y.titles[c("absolute.boxplots", "relative.boxplots") %in% devnplots]
    }
    if (length(y.titles) != length(devnplots))
      stop("y.titles does not have a title for each plot in deviations.plot")
    names(y.titles) <- devnplots
    
    ggfacet <- list()
    #Set up facet if have any
    if (all(facet.x != ".") | all(facet.y != "."))
    {
      facet.form <- facet.char2formula(facet.x,facet.y)
      if (is.null(labeller))
        ggfacet <- list(facet_grid(facet.form))
      else
        ggfacet <- list(facet_grid(facet.form, labeller = labeller))
      ggfacet <- c(ggfacet, 
                   list(theme(strip.text = element_text(size=strip.text.size, face="bold"))))
    }
    
    plotDevnBox <- function(dat, y.ytitle, ggplotFuncs = NULL, printPlot = TRUE)
    {
      plt <- ggplot(data = dat, aes_string(x = x.factor, y = "deviations"), ...) +
        theme_bw() + 
        theme(panel.grid.major = element_line(colour = "grey60", size = 0.5), 
              panel.grid.minor = element_line(colour = "grey80", size = 0.5),
              axis.title = element_text(face="bold")) +
        geom_boxplot() + 
        geom_hline(yintercept=0, colour="blue") +
        ylab(y.ytitle) + xlab(x.title) + 
        ggfacet
      
      if (!is.null(ggplotFuncs))
      {
        for(f in ggplotFuncs)
          plt <- plt + f
      }
      
      if (printPlot)
        print(plt)
      
      return(plt)
    }
    
    if ("absolute.boxplots" %in% devnplots)
    {
      dat$deviations <- dat[[observed]] - dat[[smoothed]]
      plts[["absolute"]] <- plotDevnBox(dat, y.ytitle = y.titles["absolute.boxplots"], 
                                        ggplotFuncs = ggplotFuncs, printPlot = printPlot)
    }
    if ("relative.boxplots" %in% devnplots)
    {
      dat$deviations <- (dat[[observed]] - dat[[smoothed]])/dat[[smoothed]]
      plts[["relative"]] <- plotDevnBox(dat, y.ytitle = y.titles["relative.boxplots"], 
                                        ggplotFuncs = ggplotFuncs, printPlot = printPlot)
    }
  }
  invisible(plts)
}

"plotSmoothsMedianDevns" <- function(data, response, response.smoothed = NULL, 
                                     individuals = "Snapshot.ID.Tag", times = "DAP", 
                                     trait.types = c("response", "AGR", "RGR"), 
                                     x.title = NULL, y.titles = NULL, labeller = NULL, 
                                     plots.by.med = NULL, plots.group.med = NULL, 
                                     facet.x.med = ".", facet.y.med = ".", 
                                     colour.values.med = NULL, shape.values.med = NULL, 
                                     alpha.med = 0.5, 
                                     propn.note.med = TRUE, propn.types.med = c(0.1, 0.5, 0.75), 
                                     ggplotFuncsMedDevn = NULL, printPlot = TRUE, ...)  
{
  inargs <- list(...)
  checkEllipsisArgs(c("plotSmoothsMedianDevns","plotProfiles"), inargs)
  
  if (is.null(response.smoothed))
    response.smoothed <- paste0("s", response)
  #Check that responses, response.smoothed, individuals and times are in data
  checkNamesInData(c(response, response.smoothed, individuals, times), data = data)
  
  options <- c("response", "AGR", "RGR", "all")
  traits <- options[unlist(lapply(trait.types, check.arg.values, options=options))]
  if ("all" %in% traits)
    traits <- c("response", "AGR", "RGR")
  
  #Checking of the arguments that control the plots layout
  checkLayoutArgs(data = data, plots.by.med, plots.group.med, facet.x.med, facet.y.med)
  plts.by <- plots.by.med
  plts.group <- plots.group.med
  
  #Check have a valid smooths.frame
  validsmoothsframe <- validSmoothsFrame(data)  
  if (is.character(validsmoothsframe))
    stop(validsmoothsframe)
  checkPlotsArgs(data, plts.by, plts.group, facet.x = facet.x.med, facet.y = facet.y.med)
  
  strip.text.size <- 10
  
  dat <- data
  
  #Form data.frame with just columns needed 
  #Create a factor Times that has the plotted values of x for its labels
  if (is.null(x.title))
    x.title <- times
  data[times] <- convertTimes2numeric(data[[times]])
  id.cols <- c(individuals, times) #, colour.column)
  times.factor <- ".Time.fac"
  dat[times.factor] <- dat[times]
  dat[times.factor] <- with(dat, eval(parse(text =times)))
  dat[times.factor] <- factor(unlist(dat[times.factor]), 
                              labels = unique(dat[times.factor])[order(unique(dat[[times.factor]])),])
  fac.group <- NULL
  if (!is.allnull(plts.group))
  { 
    dat$SmoothParams <- fac.mixcombine(dat, plts.group, smooth.cols = smooth.cols)
    fac.group <- "SmoothParams"
    id.cols <- c(fac.group, id.cols)
  }
  
  #Determine xfacet and fac.by
  xfacet <- facet.x.med
  if (!is.allnull(plts.by))
  {
    dat$fac.by <- fac.mixcombine(dat, plts.by, smooth.cols = smooth.cols)
    id.cols <- c("fac.by", id.cols)
  }
  if (all(xfacet != "."))
    id.cols <- c(id.cols, fac.getinFormula(xfacet))
  if (all(facet.y.med != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.y.med))
  
  if (!all(id.cols %in% names(dat)))
    stop(paste("Do not have the following required columns in data: ", 
               paste(id.cols[!(id.cols %in% names(dat))],collapse=", "), "\n", sep=""))
  
  #Set up facet
  ggfacet <- list()
  facet.cols <- NULL
  if (all(xfacet != ".") | all(facet.y.med != "."))
  {
    facet.form <- facet.char2formula(xfacet, facet.y.med)
    if (is.null(labeller))
      ggfacet <- list(facet_grid(facet.form, scales = "free_x"))
    else
      ggfacet <- list(facet_grid(facet.form, scales = "free_x", labeller = labeller))
    facet.cols <- c(xfacet, facet.y.med)
    facet.cols <- facet.cols[facet.cols != "."]
  }
  
  #Form raw and smoothed trait names
  addRates <- function(traits, response, sep = ".")
  {
    unlist(lapply(traits, 
                  function(trait, response)
                  {
                    if (!("response" %in% trait))
                      response <- paste(response, trait, sep = sep)
                    return(response)
                  }, response = response))
  }
  kresp <- addRates(traits, response = response)
  kresp.sm <- addRates(traits, response = response.smoothed)
  kresp.devn <- paste(kresp, "devn", sep = ".")
  names(kresp.sm) <- kresp
  names(kresp.devn) <- kresp
  if (is.null(y.titles))
  {
    y.titles <- paste("median", addRates(traits, response = response, sep = " "), "deviations")
    names(y.titles) <- kresp
  } else
  {
    if (length(y.titles) != length(kresp))
      stop("y.titles should be the same length as trait.types")
    else
      names(y.titles) <- kresp
  }
  
  #Calculate the deviations
  dat[kresp.devn] <- dat[kresp] - dat[kresp.sm]
  
  #Calculate the median deviations
  split.facs <- c(if (!is.allnull(plts.by)) "fac.by", 
                  facet.cols, fac.group, times.factor) #NULL objects will be ignored
  dat <- dat[c(split.facs, setdiff(names(dat), split.facs))]
  dat <- dat[do.call(order, dat),]
  tmp <- dat
  tmp$split.fac <- dae::fac.combine(as.list(dat[split.facs]), combine.levels = TRUE)
  dat.split <- split(tmp, f = tmp$split.fac)
  med.devn.dat <- lapply(dat.split, 
                         function(data, kresp.devn)
                         { 
                           med <- unlist(lapply(data[kresp.devn], median, na.rm = TRUE))
                           krow <- cbind(data[1,split.facs], rbind(med))
                           return(krow)
                         },
                         kresp.devn = kresp.devn)
  med.devn.dat <- as.data.frame(do.call(rbind, med.devn.dat))
  med.devn.dat[times] <- dae::as.numfac(unlist(med.devn.dat[times.factor]))
  #Remove the times.factor
  med.devn.dat <- med.devn.dat[,-match(times.factor, names(med.devn.dat))]
  
  #Calculate the median responses
  if (!is.null(propn.types.med))
  {
    if (length(propn.types.med) != length(kresp))
      stop("Length of propn.types.med is not the same as the number of trait.types")
    names(propn.types.med) <- kresp
    med.resp.dat <- lapply(dat.split, 
                           function(data, kresp)
                           { 
                             med <- unlist(lapply(data[kresp], median, na.rm = TRUE))
                             krow <- cbind(data[1,split.facs], rbind(med))
                             return(krow)
                           },
                           kresp = kresp)
    med.resp.dat <- as.data.frame(do.call(rbind, med.resp.dat))
    med.resp.dat <- rbind(med.resp.dat, med.resp.dat)
    med.resp.dat <- cbind(sign = rep(c(1,-1), each = nrow(med.resp.dat)/2),
                          med.resp.dat)
    med.resp.dat[kresp] <- as.data.frame(mapply(function(var, propn)
    {
      var <- propn * rep(c(1, -1), each = length(var)/2) * var
    },
    med.resp.dat[kresp], propn.types.med))
    med.resp.dat[times] <- dae::as.numfac(unlist(med.resp.dat[times.factor]))
    #Remove the times.factor
    med.resp.dat <- med.resp.dat[,-match(times.factor, names(med.resp.dat))]
  }
  
  #Plot the median deviations for each trait
  if (is.null(shape.values.med))
    shape.values.med <- c(21:24,7,9,10,11,3,4)
  plts <- list()
  for (k in kresp)
  {
    plts[[k]] <- list()
    if (is.allnull(plts.by))
      levs.by <- "all"
    else
      levs.by <- levels(med.devn.dat$fac.by)
    for (p in levs.by)
    { 
      if (is.allnull(plts.by))
        tmp <- med.devn.dat
      else
        tmp <- subset(med.devn.dat, med.devn.dat$fac.by == p)
      if ("Method" %in% names(tmp))
        tmp$Method <- with(tmp, dae::fac.recast(Method, 
                                                newlevels = substring(levels(Method),1,3)))
      plts[[k]][[p]] <- ggplot(tmp, aes_string(x = times, kresp.devn[k]), ...) +
        ggfacet +
        geom_hline(yintercept=0, linetype="solid", size=0.5, colour = "maroon", alpha=0.7) +
        scale_x_continuous(breaks = seq(min(tmp[[times]], na.rm = TRUE),
                                        max(tmp[[times]], na.rm = TRUE), by = 4)) + 
        xlab(x.title) + ylab(y.titles[k]) + theme_bw() +
        theme(strip.text = element_text(size=strip.text.size, face="bold"),
              axis.title = element_text(face="bold"),
              panel.grid.major = element_line(colour = "grey60", size = 0.5), 
              panel.grid.minor = element_line(colour = "grey80", size = 0.5))
      
      if (is.null(fac.group))
      { 
        if (is.allnull(colour.values.med))
        { 
          plts[[k]][[p]] <- plts[[k]][[p]] + geom_line (size=0.4, alpha=alpha.med)
          if (is.allnull(shape.values.med))
            plts[[k]][[p]] <- plts[[k]][[p]] + geom_point(alpha=alpha.med, size=1.5)
          else
            plts[[k]][[p]] <- plts[[k]][[p]] + 
              geom_point(shape=shape.values.med[1], alpha=alpha.med, size=1.5)
        }
        else
        { 
          plts[[k]][[p]] <- plts[[k]][[p]] + 
            geom_line (colour=colour.values.med[1], size=0.4, alpha=alpha.med)
          if (is.allnull(shape.values.med))
            plts[[k]][[p]] <- plts[[k]][[p]] + 
              geom_point(colour=colour.values.med[1], 
                         fill=colour.values.med[1], alpha=alpha.med, size=1.5)
          else
            plts[[k]][[p]] <- plts[[k]][[p]] + 
              geom_point(colour=colour.values.med[1], shape=shape.values.med[1], 
                         fill=colour.values.med[1], alpha=alpha.med, size=1.5)
        }
      }
      else
        plts[[k]][[p]] <- plts[[k]][[p]] + 
        geom_line (aes_string(colour=fac.group), size=0.4, alpha=alpha.med) +
        geom_point(aes_string(colour=fac.group, shape=fac.group, fill=fac.group), 
                   alpha=alpha.med, size=1.5)
      
      if (!(is.null(colour.values.med)))
        plts[[k]][[p]] <- plts[[k]][[p]] + scale_colour_manual(values = colour.values.med)
      
      if (!(is.null(shape.values.med)))
        plts[[k]][[p]] <- plts[[k]][[p]] + scale_shape_manual(values = shape.values.med)
      else
        plts[[k]][[p]] <- plts[[k]][[p]] + scale_shape_manual(values = c(21:24,7,9,10,11))
      
      if (!is.allnull(plts.group) && fac.group == "SmoothParams")
        plts[[k]][[p]] <- plts[[k]][[p]] + guides(shape=guide_legend(title = "Smoothing\nparameters"),
                                                  fill=guide_legend(title = "Smoothing\nparameters"),
                                                  colour=guide_legend(title = "Smoothing\nparameters"))
      
      if (!is.allnull(plts.by))
        plts[[k]][[p]] <- plts[[k]][[p]] + ggtitle(paste("Plot for",p))
      
      #Plot an envelope of the response median
      if (!is.null(propn.types.med))
      {
        if (is.allnull(plts.by))
          med.resp.tmp <- med.resp.dat
        else
          med.resp.tmp <- subset(med.resp.dat, med.resp.dat$fac.by == p)
        #Construct message to be plotted
        if (propn.note.med)
        {
          xmin <- min(med.resp.tmp[times], na.rm = TRUE)
          xrange <- max(med.resp.tmp[times], na.rm = TRUE) - xmin #not used at present
          ymin <- min(med.resp.tmp[k], na.rm = TRUE)
          envel <- data.frame(rep(xmin, 2),
                              c(ymin+2.75, ymin))
          names(envel) <- c(times, kresp.devn[k])
          ncol.x <- 1
          if (all(facet.x.med != "."))
          {
            lastlev <- levels(unlist(med.resp.tmp[facet.x.med]))
            ncol.x <- length(lastlev)
            lastlev <- factor(lastlev[length(lastlev)], levels = lastlev)
            envel[facet.x.med] <- rep(lastlev, 2)
          }
          if (all(facet.y.med != "."))
          {
            lastlevs <- lapply(facet.y.med, 
                               function(fac, med.resp.tmp) 
                               {
                                 lastlev <- levels(med.resp.tmp[[fac]])
                                 lastlev <- rep(factor(lastlev[length(lastlev)], levels = lastlev), 2)
                               }, med.resp.tmp = med.resp.tmp)
            names(lastlevs) <- facet.y.med
            envel <- cbind(envel, lastlevs)
          }
          if (ncol.x >3)
            envel$lab <- c("Envelope:", 
                           paste(propn.types.med[k],"of response median"))
          else
          {
            envel <- envel[2,]
            envel$lab <- paste("Envelope:", propn.types.med[k],"of response median")
          }
        }
        
        #Plot the envelope
        plts[[k]][[p]] <- plts[[k]][[p]] + geom_line(data = med.resp.tmp, aes_string(y=k, group="sign"), 
                                                     linetype="dashed")
        if (propn.note.med)
          plts[[k]][[p]] <- plts[[k]][[p]] + 
          geom_text(data = envel, 
                    mapping = aes_string(x = times, y = kresp.devn[k], 
                                         label = "lab"), 
                    hjust = 0, vjust=-Inf, 
                    fontface = "plain", size = 3)
      }
      
      if (!is.null(ggplotFuncsMedDevn))
      {
        for(f in ggplotFuncsMedDevn)
          plts[[k]][[p]] <- plts[[k]][[p]] + f
      }
      if (printPlot)
        print(plts[[k]][[p]])
    }
  }
  invisible(list(plots = plts, med.devn.dat = med.devn.dat))
}

"plotSmoothsComparison" <- function(data, response, response.smoothed = NULL,
                                    individuals = "Snapshot.ID.Tag", times = "DAP", 
                                    trait.types = c("response", "AGR", "RGR"), 
                                    x.title = NULL, y.titles = NULL, labeller = NULL, 
                                    which.plots = "profiles", printPlot = TRUE, 
                                    plots.include.raw = FALSE, 
                                    plots.by.pf = NULL, facet.x.pf = ".", facet.y.pf = ".", 
                                    colour.pf = "black", colour.column.pf = NULL, 
                                    colour.values.pf = NULL, alpha.pf = 0.3, 
                                    addMediansWhiskers.pf = TRUE,
                                    ggplotFuncsProfile = NULL, ggplotFuncsDevnBoxes = NULL, 
                                    ...)
{
  inargs <- list(...)
  checkEllipsisArgs(c("plotSmoothsComparison","plotProfiles"), inargs)
  
  #Find out if any plotProfiles arguments that probeSmooths handles have been supplied in '...'
  if (length(inargs))
  {
    usedProfile.args <- c("data","times","response","individuals","facet.x","facet.y",
                          "labeller","colour","colour.column","colour.values",
                          "alpha","y.title","printPlot","ggplotFuncs", "addMediansWhiskers")
    doubleargs <- intersect(names(inargs), usedProfile.args)
    if (length(doubleargs))
      stop("the  'plotProfiles' arguments ",paste0(doubleargs, collapse = ", "), 
           " conflict with 'traitSmooth' arguments")
    #extract any remaining plotProfiles arguments from inargs
    pltLongi.args <- setdiff(formalArgs(plotProfiles), c(usedProfile.args))
    pltLongi.args <- names(inargs)[names(inargs) %in% pltLongi.args]
    if (length(pltLongi.args))
      pltLongi.args <- inargs[pltLongi.args]
    else
      pltLongi.args <- NULL
  } else
    pltLongi.args <- NULL
  
  options <- c("response", "AGR", "RGR", "all")
  traits <- options[unlist(lapply(trait.types, check.arg.values, options=options))]
  if ("all" %in% traits)
    traits <- c("response", "AGR", "RGR")
  
  #Checking of the arguments that control the plots layout
  checkLayoutArgs(data = data, plots.by.pf, plts.group = NULL, facet.x.pf, facet.y.pf)
  plts.by <- plots.by.pf
  
  #Check have a valid smooths.frame
  validsmoothsframe <- validSmoothsFrame(data)  
  if (is.character(validsmoothsframe))
    stop(validsmoothsframe)
  checkPlotsArgs(data, plts.by, facet.x = facet.x.pf, facet.y = facet.y.pf)
  
  options <- c("none", "profiles", "absolute.boxplots", "relative.boxplots", "medians.deviations")
  plots <- options[unlist(lapply(which.plots, check.arg.values, options=options))]
  if ("none" %in% plots & length(plots) > 1)
    plots <- "none"
  
  if (is.null(x.title))
    x.title <- times
  if (is.null(response.smoothed))
    response.smoothed <- paste0("s", response)
  response.smooth <- response.smoothed
  #Check that responses, response.smoothed, individuals and times are in data
  checkNamesInData(c(response, response.smoothed, individuals, times), data = data)
  
  addRates <- function(traits, response, sep = ".")
  {
    unlist(lapply(traits, 
                  function(trait, response)
                  {
                    if (!("response" %in% trait))
                      response <- paste(response, trait, sep = sep)
                    return(response)
                  }, response = response))
  }
  kresp <- addRates(traits, response = response)
  if (!all(kresp %in% names(data)))
    stop("The following traits are not in the smooths.frame: ",
         paste0(kresp[!(kresp %in% names(data))], collapse = ", "), 
         "; perhaps, trait.types needs to be set differently")
  kresp.sm <- addRates(traits, response = response.smoothed)
  names(kresp.sm) <- kresp
  
  if (is.null(y.titles))
  {
    y.titles <- addRates(traits, response = response, sep = " ")
    names(y.titles) <- kresp
  } else
  {
    if (length(y.titles) != length(kresp))
      stop("y.titles should be the same length as trait.types")
    else
      names(y.titles) <- kresp
  }
  
  data[times] <- convertTimes2numeric(data[[times]])
  times.factor <- ".Time.fac"
  data[times.factor] <- data[times]
  data[times.factor] <- with(data, eval(parse(text =times)))
  data[times.factor] <- factor(unlist(data[times.factor]), 
                               labels = unique(data[times.factor])[order(unique(data[[times.factor]])),])
  
  #Set up the facets  
  if (all(facet.x.pf == "."))
    xfacet = facet.x.pf
  else
  {
    data$SmoothParams <- fac.mixcombine(data, facet.x.pf, smooth.cols = smooth.cols)
    xfacet <- "SmoothParams"
  }
  
  #Do the plots
  x.axis <- list(scale_x_continuous(breaks = seq(min(data[[times]], na.rm = TRUE),
                                                 max(data[[times]], na.rm = TRUE), 
                                                 by = 4)),
                 theme(axis.text.x = element_text(size = 7.5)))
  plts <- list()
  for (k in kresp)
  {
    plts[[k]] <- list()
    plts[[k]][["deviations"]] <- plts[[k]][["profiles"]] <- list()
    if (all(facet.x.pf == ".") && plots.include.raw &&  ("profiles" %in% plots))
    { 
      plts[[k]][["profiles"]][["Unsmoothed"]] <- 
        do.call(plotProfiles, 
                c(list(data = data, times = times, response = k, 
                       individuals = individuals, 
                       facet.x=xfacet, facet.y=facet.y.pf, 
                       labeller = labeller, 
                       colour = colour.pf, 
                       colour.column = colour.column.pf, 
                       colour.values = colour.values.pf, 
                       alpha = alpha.pf, 
                       title="Plot of unsmoothed response", 
                       x.title = x.title, y.title = y.titles[k], 
                       addMediansWhiskers = addMediansWhiskers.pf, 
                       printPlot=FALSE, 
                       ggplotFuncs = c(ggplotFuncsProfile, x.axis)), 
                  pltLongi.args))
      if (printPlot)
        print(plts[[k]][["profiles"]][["Unsmoothed"]])
    }
    if (is.allnull(plts.by))
      levs.by <- "all"
    else
    {
      data$plots.by.pf <- fac.mixcombine(data, plts.by, smooth.cols = smooth.cols)
      levs.by <- levels(data$plots.by.pf)
    }
    #Loop over plots.by.pf
    for (by in levs.by)
    {
      if ("profiles" %in% plots)
      { 
        if (is.allnull(plts.by))
        { 
          title <- NULL
          tmp1 <- data
        } else
        {
          title <- paste0("Plot for ", by)
          tmp1 <- data[data$plots.by.pf==by,]
        }
        if (!(xfacet == ".") && plots.include.raw)
        {
          tmp2 <- tmp1
          tmp2[kresp.sm[k]] <- tmp2[k]
          tmp2["SmoothParams"] <- "Raw"
          levs <- c("Raw", levels(factor(tmp1[["SmoothParams"]])))
          tmp1 <- rbind(tmp2,tmp1)
          tmp1["SmoothParams"] <- factor(tmp1[["SmoothParams"]], levels = levs) 
        }
        plts[[k]][["profiles"]][[by]] <- 
          do.call(plotProfiles, 
                  c(list(data = tmp1, times = times, 
                         response = kresp.sm[k], 
                         individuals = individuals, 
                         facet.x=xfacet, facet.y=facet.y.pf, 
                         labeller = labeller, 
                         colour = colour.pf, 
                         colour.column = colour.column.pf, 
                         colour.values = colour.values.pf, 
                         alpha = alpha.pf, 
                         title = title, 
                         x.title = x.title, y.title = y.titles[k], 
                         addMediansWhiskers = addMediansWhiskers.pf, 
                         printPlot=FALSE, 
                         ggplotFuncs = c(ggplotFuncsProfile, x.axis)), 
                    pltLongi.args))
        if (printPlot)
          print(plts[[k]][["profiles"]][[by]])
      }
      if (any(c("absolute.boxplots", "relative.boxplots") %in% plots))
      {
        y.titles.devn <- c(paste("Absolute", k, "deviations", sep = " "),
                           paste("Relative", k, "deviations", sep = " "))
        names(y.titles.devn) <- c("absolute.boxplots", "relative.boxplots")
        y.titles.devn <- y.titles.devn[c("absolute.boxplots", "relative.boxplots") %in% plots]
        
        #Plot deviation plots for current plots.by.pf
        if (is.allnull(plts.by))
          tmp1 <- data
        else
        { 
          tmp1 <- data[data$plots.by.pf==by,]
          ggplotFuncsDevnBoxes <- c(ggplotFuncsDevnBoxes, list(ggtitle(paste0("Plot for ", by))))
        }
        plt <- plotDeviationsBoxes(data = tmp1, x.factor = times.factor, 
                                   observed = k, smoothed = kresp.sm[k], 
                                   deviations.plots = plots, 
                                   x.title = x.title, y.titles = y.titles.devn, 
                                   facet.x=xfacet, facet.y=facet.y.pf, 
                                   df = degfree, ggplotFuncs = ggplotFuncsDevnBoxes,
                                   printPlot = printPlot)
        plts[[k]][["deviations"]][["absolute"]][[by]] <- plt[["absolute"]]
        plts[[k]][["deviations"]][["relative"]][[by]] <- plt[["relative"]]
      }
    }    
  }
  invisible(plts)
}

#Function to fit splines, including possible boundary correction from Huang (2001)
ncsSpline <- function(vars, correctBoundaries = FALSE, 
                      df, lambda, cv = FALSE,  ...)
{
  if (ncol(vars) != 2)
    stop("Must supply a two-column matrix or data.frame")
  if (!correctBoundaries)
  {
    if (missing(df))
    {
      if (missing(lambda))
        fity <- smooth.spline(vars, all.knots=TRUE, ...)
      else
        fity <- smooth.spline(vars, all.knots=TRUE, lambda = lambda, ...)
      
    } else
    {
      if (missing(lambda))
        fity <- smooth.spline(vars, all.knots=TRUE, df=df, ...)
      else
        stop("Only one of df and lambda can be specified")
    }
    fit.spline <- list(x = fity$x, 
                       y = fity$y, 
                       lev = fity$lev,
                       lambda = fity$lambda,
                       df = fity$df,
                       uncorrected.fit = fity)  
  } else
  {
    nval <- nrow(vars)
    W <- matrix(NA, nrow = nval, ncol = 4)
    W2 <- matrix(NA, nrow = nval, ncol = 4)
    W3 <- matrix(NA, nrow = nval, ncol = 4)
    x <- vars[,1]
    y <- vars[,2]
    
    # construct the four polynomials
    W[,1] <- x^2/2-x^4/4+x^5/10
    W[,2] <- x^3/6-x^4/6+x^5/20
    W[,3] <- x^4/4-x^5/10
    W[,4] <- -x^4/12+x^5/20
    if (missing(df))
    {
      if (!missing(lambda))
        stop("lambda must not be set for correctBoundaries = TRUE")
      lam <- vector(mode = "numeric", length = nval)
      rss <- vector(mode = "numeric", length = nval)
      tr <- vector(mode = "numeric", length = nval)
      gcv <- vector(mode = "numeric", length = nval)
      # use GCV to search for best lambda
      for(i in 1:nval) 
      {
        u <- -7+7.0*(i-1)/(nval-1)
        lam[i] <- 10^u
        slam <- lam[i]
        # get regular fit for y and four polynomials
        fity <- smooth.spline(x,y,all.knots=TRUE,spar=slam, ...)
        W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,spar=slam, ...)$y
        W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,spar=slam, ...)$y
        W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,spar=slam, ...)$y
        W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,spar=slam, ...)$y
        #Apply boundary correction to the fitted spline
        W2 <- W-W2
        h <- solve(t(W2)%*%W2,t(W2)%*%(y-fity$y))
        # get the newfit, rss and calculate the trace of new S
        newfit <- fity$y+W2%*%h
        rss[i] <- sum((newfit-y)^2)
        W3[,1] <- smooth.spline(x,W2[,1],all.knots=TRUE,spar=slam, ...)$y
        W3[,2] <- smooth.spline(x,W2[,2],all.knots=TRUE,spar=slam, ...)$y
        W3[,3] <- smooth.spline(x,W2[,3],all.knots=TRUE,spar=slam, ...)$y
        W3[,4] <- smooth.spline(x,W2[,4],all.knots=TRUE,spar=slam, ...)$y
        W3 <- W2-W3
        K <- solve(t(W2)%*%W2,t(W3))
        tr[i] <- sum(fity$lev)+sum(diag(K%*%W2))
        gcv[i] <- nval*rss[i]/(nval-tr[i])^2
      }
      # get the optimal lambda and apply to y and the four polynomials
      lam.opt<- lam[order(gcv)[1]]
      slam <- lam.opt
      fity <- smooth.spline(x,y,all.knots=TRUE,spar=slam)
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,spar=slam, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,spar=slam, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,spar=slam, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,spar=slam, ...)$y
    } else #df is specified
    {
      # get regular fit for y and four polynomials for specified df
      fity <- smooth.spline(x,y,all.knots=TRUE,df=df, ...)
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,df=df, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,df=df, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,df=df, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,df=df, ...)$y
      # calculate the trace of new S
      W3[,1] <- smooth.spline(x,W2[,1],all.knots=TRUE,df=df, ...)$y
      W3[,2] <- smooth.spline(x,W2[,2],all.knots=TRUE,df=df, ...)$y
      W3[,3] <- smooth.spline(x,W2[,3],all.knots=TRUE,df=df, ...)$y
      W3[,4] <- smooth.spline(x,W2[,4],all.knots=TRUE,df=df, ...)$y
      W3 <- W2-W3
      K <- solve(t(W2)%*%W2,t(W3))
    }
    #Apply boundary correction to the fitted spline
    W2 <- W-W2
    h <- solve(t(W2)%*%W2,t(W2)%*%(y-fity$y))
    newy <- as.vector(fity$y+W2%*%h)
    # get the newfit leverage values of new S
    lev <- as.vector(fity$lev+diag(W2%*%K))
    tr <- sum(lev)
    rss <- sum((y - newy)^2)
    lambda <- nval*rss/(nval-tr)^2
    # Set up list to return
    fit.spline <- list(x = fity$x, 
                       y = newy, 
                       lev = lev,
                       lambda = lambda,
                       df = fity$df,
                       uncorrected.fit = fity)  
    if (!missing(df))
      fit.spline$df <- df
  }
  class(fit.spline) <- "ncsSpline"
  return(fit.spline)
}

predict.ncsSpline <- function(object, x, correctBoundaries = FALSE, 
                              df, cv = FALSE,  ...)
{
  if (class(object) != "ncsSpline")
    stop("Must supply a an object of class ncsSpline")
  fit <- predict(object$uncorrected.fit, x = x)
  
  #Correct boundaries
  if (correctBoundaries)
  {
    nval <- length(x)
    W <- matrix(NA, nrow = nval, ncol = 4)
    W2 <- matrix(NA, nrow = nval, ncol = 4)
    vars <- data.frame(x = object$uncorrected.fit$x, 
                       yin = object$uncorrected.fit$yin)
    vars <- merge(as.data.frame(fit), vars, all.x = TRUE)
    vars$yin[is.na(vars$yin)] <- vars$y[is.na(vars$yin)]
    vars <- vars[order(vars$x), ]
    x <- vars$x
    
    # construct the four polynomials
    W[,1] <- x^2/2-x^4/4+x^5/10
    W[,2] <- x^3/6-x^4/6+x^5/20
    W[,3] <- x^4/4-x^5/10
    W[,4] <- -x^4/12+x^5/20
    if (missing(df))
    {
      slam <- object$lambda
      # get regular fit for y and four polynomials
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,spar=slam, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,spar=slam, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,spar=slam, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,spar=slam, ...)$y
    } else #df is specified
    {
      # get regular fit for y and four polynomials for specified df
      W2[,1] <- smooth.spline(x,W[,1],all.knots=TRUE,df=df, ...)$y
      W2[,2] <- smooth.spline(x,W[,2],all.knots=TRUE,df=df, ...)$y
      W2[,3] <- smooth.spline(x,W[,3],all.knots=TRUE,df=df, ...)$y
      W2[,4] <- smooth.spline(x,W[,4],all.knots=TRUE,df=df, ...)$y
    }
    #Apply boundary correction to the fitted spline
    W2 <- W-W2
    h <- solve(t(W2)%*%W2,t(W2)%*%(vars$yin-vars$y))
    fit <- list(x = vars$x, 
                y = as.vector(vars$y+W2%*%h))  
  }
  return(fit)
}

#Functions to fit P-splines using JOPS
pSpline <- function(vars, npspline.segments, lambda = NULL, ...)
{
  if (ncol(vars) != 2)
    stop("Must supply a two-column matrix or data.frame")
  fity <- JOPS::psNormal(x = vars[[1]], y = vars[[2]], nseg = npspline.segments, lambda = lambda, 
                         xgrid = vars[[1]])
  fit.spline <- list(x = fity$x, 
                     y = as.vector(fity$muhat), 
                     lev = NULL,
                     lambda = fity$lambda,
                     df = fity$effdim,
                     npspline.segments = npspline.segments, 
                     uncorrected.fit = fity)  
  class(fit.spline) <- "PSpline"
  return(fit.spline)
}

predict.pSpline <- function(object, x, npspline.segments, deriv = 0)
{
  fit.obj <- object$uncorrected.fit
  bdeg <- fit.obj$bdeg
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (class(object) != "PSpline")
    stop("Must supply a an object of class PSpline")
  if (deriv == 0)
  {
    preds <- predict(fit.obj, x = x, type = "mu")
    fit <- list(x = x, y = preds)
  } else #obtain a derivative of order deriv based on Eqn 2.15 from Eilers and Marx (2021) JOPS
  {
    if ((bdeg - deriv) <= 0)
      stop("The degree of the spline (3) is insufficient to compute a derivative of order ", deriv)
    alpha <- as.vector(fit.obj$pcoef)
    alphaDeriv <- diff(alpha, differences = deriv) / (((xmax - xmin)/fit.obj$nseg)^deriv)
    Bderiv <- JOPS::bbase(x, xl = xmin, xr = xmax, nseg = fit.obj$nseg, bdeg = (bdeg - deriv))
    fit <- as.vector(Bderiv %*% alphaDeriv)
    fit <- list(x = x, y = fit)
  }
  return(fit)
}

#Function to fit a spline using smooth.spline or JOPS
"smoothSpline" <- function(data, response, response.smoothed = NULL, x, 
                           smoothing.method = "direct", 
                           spline.type = "NCSS",  df=NULL, lambda = NULL, 
                           npspline.segments = NULL, correctBoundaries = FALSE, 
                           rates = NULL, suffices.rates = NULL, 
                           extra.derivs = NULL, suffices.extra.derivs=NULL, 
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
  
  #check input arguments
  impArgs <- match.call()
  if ("na.rm" %in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  if ("smoothing.scale" %in% names(impArgs))
    stop("smoothing.scale has been deprecated; use smoothing.method")
  
  #Check that required cols are in data
  checkNamesInData(c(response, x), data = data)
  
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
  
  options <- c("AGR", "PGR", "RGR")
  if (is.allnull(rates))
    grates <- NULL
  else
    grates <- options[unlist(lapply(rates, check.arg.values, options=options))]
  if (correctBoundaries && !is.allnull(grates))
    stop("Unable to correctBoundaries when rates is not NULL")
  if (!is.allnull(grates) && "PGR" %in% grates)
    stop("PGR is not available when rates are based on derivatives")
  if (!is.null(suffices.rates) & length(grates) != length(suffices.rates))
    stop("The length of of rates and suffices.rates should be equal")
  if (!is.allnull(grates) && !is.null(extra.derivs) && (1 %in% extra.derivs))
    stop("when rates is not NULL, 1 should not be included in extra.derivs")
  if (is.null(suffices.rates))
    suffices.rates <- grates
  names(suffices.rates) <- grates
  
  if (!is.null(extra.derivs) & !is.null(suffices.extra.derivs))
    if (length(extra.derivs) != length(extra.derivs))
      stop("The number of names supplied must equal the number of derivatives specified")

  #Determine what is required from spline fitting
  if (!is.allnull(grates))
  {
    derivs = 1
    if (smethod == "direct")
    {
      if (!("AGR" %in% grates))
        suffices.derivs <- "_tmp"
      else
      {
        if (is.null(suffices.rates))
          suffices.derivs <- "AGR"
        else
          suffices.derivs <- suffices.rates["AGR"]
      }
      names(suffices.derivs) <- "AGR"
      if ("RGR" %in% grates)
      { 
        extra.rate <- "RGR"
        names(extra.rate) <- suffices.rates["RGR"]
      }
      else
        extra.rate <- NULL
    } else
    {
      if (!("RGR" %in% grates))
        suffices.derivs <- "_tmp"
      else
      {
        if (is.null(suffices.rates))
          suffices.derivs <- "RGR"
        else
          suffices.derivs <- suffices.rates["RGR"]
      }
      names(suffices.derivs) <- "RGR"
      if ("AGR" %in% grates)
      { 
        extra.rate <- "AGR"
        names(extra.rate) <- suffices.rates["AGR"]
      }
      else
        extra.rate <- NULL
    }
  } else
  {  
    derivs = NULL
    suffices.derivs <- NULL
    extra.rate <- NULL
  }
  
  if (!is.null(extra.derivs))
  {  
    derivs <- unique(c(derivs, extra.derivs))
    suffices.derivs <- c(suffices.derivs, suffices.extra.derivs)
    names(derivs) <- suffices.derivs
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
    response.smoothed <- paste0("s", response)
  fit.names <- c(x, response.smoothed)
  if (!is.allnull(derivs))
  {
    if (is.allnull(suffices.derivs))
      fit.names <- c(fit.names, paste(response.smoothed,".dv",derivs,sep=""))
    else
      fit.names <- c(fit.names, paste(response.smoothed, suffices.derivs, sep="."))
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
    if (!is.allnull(derivs))
    {
      if (is.allnull(suffices.derivs))
        fit.names <- c(fit.names, paste(response.smoothed,".dv",derivs,sep=""))
      else
        fit.names <- c(fit.names, paste(response.smoothed, suffices.derivs, sep="."))
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
        fit <- predict.pSpline(fit.spline, x = x.pred, 
                               npspline.segments = fit.spline$uncorrected.fit$npspline.segments)
    }
    rsmooth <- response.smoothed
    names(fit) <- c(x, rsmooth)
    #backtransform if transformed
    if (smethod == "logarithmic")
      fit[[rsmooth]] <- exp(fit[[rsmooth]])
    
    #get derivatives if required
    if (!correctBoundaries & !is.null(derivs))
    {
      for (d in derivs)
      {
        if (is.null(suffices.derivs))
          rsmooth.dv <- paste0(response.smoothed,".dv",d)
        else
        { 
          k <- ifelse(length(derivs) == 1, 1, match(d, derivs))
          rsmooth.dv <- paste(response.smoothed, suffices.derivs[k], sep=".")
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
        if (is.null(suffices.derivs))
          rsmooth.dv <- paste0(response.smoothed,".dv",1)
        else
          rsmooth.dv <- paste(response.smoothed, suffices.derivs["AGR"], sep=".")
        #get the extra derivative
        if (!(rsmooth.dv %in% names(fit)))
          stop("First derivative not available to calculate RGR")
        fit[[paste(rsmooth,names(extra.rate),sep=".")]] <- fit[[rsmooth.dv]]/fit[[rsmooth]]
      }
      #Add AGR if required
      if (!is.null(extra.rate) && extra.rate == "AGR")
      { 
        #Check have the required computed derivative 
        if (is.null(suffices.derivs))
          rsmooth.dv <- paste0(response.smoothed,".dv",1)
        else
          rsmooth.dv <- paste(response.smoothed, suffices.derivs["RGR"], sep=".")
        #get the extra derivative
        if (!(rsmooth.dv %in% names(fit)))
          stop("First derivative not available to calculate AGR")
        fit[[paste(rsmooth,names(extra.rate),sep=".")]] <- fit[[rsmooth.dv]]*fit[[rsmooth]]
      }
    }
    fit <- as.data.frame(fit)
    #Remove temporary rate, if there are any
    if (length(grep("._tmp", names(fit), fixed = TRUE)))
      fit <- fit[, -grep("._tmp", names(fit), fixed = TRUE)]
    
    
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

#Function to perform the smooths specified by a set of smoothing-parameter schemes
smoothSchemes <- function(tmp, spar.schemes, 
                          response, response.smooth, times, individuals, 
                          traits, get.rates, ratemeth.opt, grates, 
                          ntimes2span = 2, nseg, correctBoundaries, 
                          na.x.action, na.y.action)
{
  
  #Form smoothed traits and grates for each scheme
  smth <- list()
  for (k in 1:nrow(spar.schemes))
  {
    scheme.lab <- paste(spar.schemes[k,], collapse = "-")
    scheme.type <- spar.schemes$Type[k]
    scheme.tpar <- spar.schemes$TunePar[k]
    scheme.tval <- spar.schemes$TuneVal[k]
    scheme.meth <- spar.schemes$Method[k]
    if (scheme.type == "NCSS")
    {
      if (scheme.tpar == "df")
      {
        scheme.df <- scheme.tval
        scheme.lambda <- NULL
      } else
      {
        scheme.df <- NULL
        scheme.lambda <- scheme.tval
      }
    } else #PS
    {
      scheme.df <- NULL
      scheme.lambda <- scheme.tval
    }
    
    rates.meth <- ratemeth.opt
    if (get.rates)
      responses.smooth <- c(response.smooth, 
                            paste(response.smooth, grates, sep="."))
    else
      rates.meth <- "none"
    
    smth[[scheme.lab]] <- byIndv4Times_SplinesGRs(data = tmp, response, response.smoothed = response.smooth, 
                                                  individuals = individuals, times = times, 
                                                  rates.method = rates.meth, which.rates = grates, 
                                                  smoothing.method = scheme.meth, 
                                                  spline.type = scheme.type, 
                                                  df= scheme.df, lambda = scheme.lambda, 
                                                  npspline.segments = nseg, 
                                                  correctBoundaries = correctBoundaries, 
                                                  na.x.action = na.x.action, 
                                                  na.y.action = na.y.action)
    smth[[scheme.lab]] <- cbind(Type = scheme.type, TunePar = scheme.tpar, TuneVal = scheme.tval, 
                                Tuning = paste(scheme.tpar, scheme.tval, sep = "-"), 
                                Method = scheme.meth, 
                                smth[[scheme.lab]])
  }
  
  #Form single data.frame 
  smth <- do.call(rbind, smth)
  cols <- names(smth)
  resps <- names(smth)[grepl(response, names(smth), fixed = TRUE)]
  cols <- setdiff(cols,resps)
  smth <- smth[c(cols,resps)]
  rownames(smth) <- NULL
  #Convert smoothing combinations to factors, paying attention to levels order 
  smth[c("Type","TunePar","TuneVal","Tuning","Method")] <- 
    mapply(function(x, sch) factor(x, levels = as.character(unique(sch))), 
           smth[c("Type","TunePar","TuneVal","Tuning","Method")], spar.schemes, SIMPLIFY = FALSE)
  
  return(smth)
}

"probeSmooths" <- function(data, response = "PSA", response.smoothed = NULL, 
                           individuals="Snapshot.ID.Tag", times = "DAP", 
                           keep.columns = NULL, 
                           trait.types = c("response", "AGR", "RGR"), 
                           get.rates = TRUE, rates.method="differences", 
                           ntimes2span = NULL, 
                           smoothing.methods = "direct", smoothing.segments = NULL, 
                           spline.types = "NCSS", df=NULL, lambdas = NULL, 
                           npspline.segments = NULL, smoothing.schemes = NULL, 
                           na.x.action="exclude", na.y.action = "trimx", 
                           external.smooths = NULL, 
                           correctBoundaries = FALSE, 
                           x.title = NULL, y.titles = NULL, labeller = NULL, 
                           which.plots = "profiles", 
                           plots.include.raw = FALSE,
                           plots.by.pf = NULL, facet.x.pf = ".", facet.y.pf = ".", 
                           colour.pf = "black", colour.column.pf = NULL, 
                           colour.values.pf = NULL, alpha.pf = 0.3, 
                           addMediansWhiskers.pf = TRUE,
                           ggplotFuncsProfile = NULL, 
                           plots.by.med = NULL, plots.group.med = NULL, 
                           facet.x.med = ".", facet.y.med = ".",
                           colour.values.med = NULL, shape.values.med = NULL, 
                           alpha.med = 0.5, 
                           propn.note.med = TRUE, propn.types.med = c(0.1, 0.5, 0.75), 
                           ggplotFuncsMedDevn = NULL, 
                           ggplotFuncsDevnBoxes = NULL, ...)
{ 
  #check input arguments
  impArgs <- match.call()
  if ("na.rm" %in% names(impArgs))
    stop("na.rm has been deprecated; use na.x.action and na.y.action")
  if ("smoothing.scales" %in% names(impArgs))
    stop("smoothing.scales has been deprecated; use smoothing.methods")
  if ("deviations.boxplots" %in% names(impArgs))
    stop("deviations.boxplots has been deprecated; use which.plots")
  inargs <- list(...)
  checkEllipsisArgs(c("probeSmooths","plotProfiles"), inargs)
  
  
  smooth.cols <- c("Type","TunePar","TuneVal","Tuning","Method")
  
  options <- c("differences","derivatives")
  ratemeth.opt <- options[check.arg.values(rates.method, options=options)]
  if (get.rates && is.null(ntimes2span))
  {
    if (ratemeth.opt == "differences")
      ntimes2span <- 2
    if (ratemeth.opt == "derivatives")
      ntimes2span <- 3
  }
  
  options <- c("none", "profiles", "absolute.boxplots", "relative.boxplots", "medians.deviations")
  plots <- options[unlist(lapply(which.plots, check.arg.values, options=options))]
  if ("none" %in% plots & length(plots) > 1)
    plots <- "none"
  if (is.null(x.title))
    x.title <- times
  data[times] <- convertTimes2numeric(data[[times]])
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
    if (!all(diff(unlist(smoothing.segments)) > 0))
      stop("the smoothing.segments are not a set of non-overlapping, successive intervals")
  }
  
  #Checking of the arguments that control the plots layout for profiles plots
  if (any(c("profiles", "absolute.boxplots", "relative.boxplots") %in% plots))
    checkLayoutArgs(data = NULL, plots.by.pf, plts.group = NULL, facet.x.pf, facet.y.pf)
  plts.by <- plots.by.pf
  
  #Checking of the arguments that control the plots layout for the medians deviations plots 
  if ("medians.deviations" %in% plots)
    checkLayoutArgs(data = NULL, plots.by.med, plots.group.med, facet.x.med, facet.y.med)
  
  plts.by.med <- plots.by.med
  plts.group.med <- plots.group.med
  smethods.opt <- c("direct", "logarithmic")
  smethods <- smethods.opt[unlist(lapply(smoothing.methods, check.arg.values, options=smethods.opt))]
  methlabs <- c("Direct", "Log")
  names(methlabs) <- smethods.opt
  methlabs <- methlabs[smethods]
  stypes.opt <- c("NCSS", "PS")
  stypes <- stypes.opt[unlist(lapply(spline.types, check.arg.values, options=stypes.opt))]
  
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
        propn.types.med <- propn.types.med[1]
        warning("get.rates is FALSE; trait.types changed to response and propn.types.med reduced to its first value")
      }
    }
    else
      grates <- c("AGR","RGR")[c("AGR","RGR") %in% traits]
  }
  
  #Form data.frame with columns needed 
  id.cols <- c(individuals, times, keep.columns, colour.column.pf)
  if (all(facet.x.pf != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.x.pf))
  if (all(facet.x.med != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.x.med))
  if (all(facet.y.pf != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.y.pf))
  if (all(facet.y.med != "."))
    id.cols <- c(id.cols, fac.getinFormula(facet.y.med))
  v <- unique(c(id.cols, response))
  v <- setdiff(v, c("Type","TunePar","TuneVal","Tuning","Method")) #remove names yet to come
  
  #Check that required cols are in data
  checkNamesInData(v, data = data)
  tmp <- data[v]
  times.diffs.in.data <- paste0(times, ".diffs") %in% names(data)
  
  #Smooth response and form growth rates
  if (is.null(response.smoothed))
    response.smooth <- paste0("s", response)
  else
    response.smooth <- response.smoothed
  responses.smooth <- response.smooth
  #Always get rates if get.rates is TRUE so that they are in the returned data
  #Do for all data even if segmented so that only the observation for the very first time is NA  
  if (get.rates) 
    tmp <- byIndv4Times_GRsDiff(data = tmp, response, 
                                individuals=individuals, 
                                times=times, avail.times.diffs = FALSE, 
                                which.rates = grates, ntimes2span = ntimes2span)
  
  #Construct combinations of smoothing parameters
  if (is.null(smoothing.schemes))
  {
    if (!is.allnull(lambdas))
    {  
      if  (!is.list(lambdas))
      {  
        warning("lambdas have been converted to a list with the single component names PS")
        lambdas <- list(PS = lambdas)
      }
      if (!all(names(lambdas) %in% stypes))
        stop("The following names for the components of lambdas are not in the specified spline.types: ",
             paste0(names(lambdas)[!(names(lambdas) %in% stypes)]))
    }
    spar.schemes.NCSS <- spar.schemes.PS <- data.frame()
    for (stype in stypes)
    {
      if (stype == "NCSS")
      {
        if (is.allnull(df) && (is.allnull(lambdas) || !("NCSS" %in% names(lambdas))))
          stop("one and only one of df and lambda must be specified for spline.type NCSS")
        if (!is.allnull(df))
        {
          spar.schemes.NCSS <-  do.call(rbind, lapply(df, 
                                                      function(kdf, TTlabs) 
                                                        cbind(TTlabs, TuneVal = kdf, 
                                                              Tuning = paste(TTlabs$TunePar, kdf, 
                                                                             sep = "-")), 
                                                      TTlabs = data.frame(Type = "NCSS",
                                                                          TunePar = "df")))
        }
        if (!is.allnull(lambdas) && ("NCSS" %in% names(lambdas))) #add lambdas if specified
        { 
          
          spar.schemes.NCSS <-  rbind(spar.schemes.NCSS, 
                                      do.call(rbind, lapply(lambdas[["NCSS"]], 
                                                            function(lambda, TTlabs) 
                                                              cbind(TTlabs, 
                                                                    TuneVal = lambda, 
                                                                    Tuning = paste(TTlabs$TunePar, lambda, 
                                                                                   sep = "-")), 
                                                            TTlabs = data.frame(Type = "NCSS",
                                                                                TunePar = "lambda"))))
        }
        #Expand to include Methods
        spar.schemes.NCSS <- do.call(rbind, lapply(smethods, 
                                                   function(smethod, spar.schemes.PS) 
                                                     cbind(spar.schemes.NCSS, Method = smethod), 
                                                   spar.schemes.PS = spar.schemes.PS))
      } else
      {
        if (stype == "PS")
        {
          if (is.allnull(lambdas) || !("PS" %in% names(lambdas)))
            stop("At least one value needs to be specfied for lambdas when spine.type is PS")
          spar.schemes.PS <-  do.call(rbind, lapply(lambdas[["PS"]], 
                                                    function(lambda, TTlabs) 
                                                      cbind(TTlabs, 
                                                            TuneVal = lambda, 
                                                            Tuning = paste(TTlabs$TunePar, lambda, 
                                                                           sep = "-")), 
                                                    TTlabs = data.frame(Type = "PS",
                                                                        TunePar = "lambda")))
          spar.schemes.PS <- do.call(rbind, lapply(smethods, 
                                                   function(smethod, spar.schemes.PS) 
                                                     cbind(spar.schemes.PS, Method = smethod), 
                                                   spar.schemes.PS = spar.schemes.PS))
        } else
          stop("unknown spline type")
      }
    }
    spar.schemes <- rbind(spar.schemes.NCSS, spar.schemes.PS)
    rownames(spar.schemes) <- NULL
    
  } else
  { 
    spar.schemes <- smoothing.schemes
    spar.schemes$Type <- stypes.opt[unlist(lapply(as.list(spar.schemes$Type), 
                                                  check.arg.values, options=stypes.opt))]
    spar.schemes$Method <- smethods.opt[unlist(lapply(as.list(spar.schemes$Method), 
                                                      check.arg.values, options=smethods.opt))]
    TPar.opt <- c("df","lambda")
    spar.schemes$TunePar <- TPar.opt[unlist(lapply(as.list(spar.schemes$TunePar), 
                                                   check.arg.values, options=TPar.opt))]
    spar.schemes$Tuning <- with(spar.schemes, paste(TunePar, TuneVal, sep = "-"))
    spar.schemes <- spar.schemes[smooth.cols]
    attr(spar.schemes, which = "nschemes") <- nrow(spar.schemes)
  }
  
  #Generate the smooths
  if (is.allnull(smoothing.segments))
    smth <- smoothSchemes(tmp = tmp, spar.schemes = spar.schemes,
                          response = response, response.smooth = response.smooth, 
                          times=times, ntimes2span = ntimes2span, 
                          individuals = individuals, traits = traits, 
                          get.rates = get.rates, ratemeth.opt = ratemeth.opt, grates = grates, 
                          nseg = npspline.segments, correctBoundaries = correctBoundaries, 
                          na.x.action = na.x.action, na.y.action = na.y.action)
  else    
  {
    knseg <- npspline.segments[1]
    smth <- data.frame()
    for (k in 1:length(smoothing.segments))
    {
      segm <- smoothing.segments[[k]]
      subdat <- tmp[(tmp[times] >= segm[1]) & (tmp[times] <= segm[2]),] 
      if (length(npspline.segments) > 1) knseg <- npspline.segments[k]
      #only get smooths when difference growth rates are required and ntimes2span is 2 
      smth <- rbind(smth, 
                    smoothSchemes(tmp = subdat, spar.schemes = spar.schemes,
                                  response = response, response.smooth = response.smooth, 
                                  times=times, ntimes2span = ntimes2span, 
                                  individuals = individuals, traits = traits, 
                                  get.rates = ifelse((ntimes2span == 2 && get.rates && 
                                                        ratemeth.opt == "differences"),
                                                     FALSE, get.rates), 
                                  ratemeth.opt = ratemeth.opt, grates = grates, 
                                  nseg = knseg, correctBoundaries = correctBoundaries, 
                                  na.x.action = na.x.action, na.y.action = na.y.action))
    }
    smth <- smth[do.call(order, smth), ]
    #get overall difference growth rates for ntimes2apn == 2
    if (ntimes2span == 2 && get.rates && ratemeth.opt == "differences") 
    {
      smth <- split(smth, fac.combine(as.list(smth[c("Type","Tuning","Method")])))
      smth <- lapply(smth, byIndv4Times_GRsDiff, responses = response.smooth, 
                     individuals = individuals, 
                     times=times, ntimes2span = ntimes2span, 
                     which.rates = grates)
      smth <- do.call(rbind, smth)
      smth <- smth[do.call(order, smth), ]
    } 
  }
  
  #Add external.smooths, if required
  if (!is.null(external.smooths))
  { 
    #Determine which smoothing-parameter columns are in external.smooths & check have required other columns 
    external.smooths <- as.data.frame(external.smooths)
    smth.cols <- intersect(smooth.cols, names(external.smooths))
    extra.vars <- setdiff(names(smth), smooth.cols)
    if (all(extra.vars %in% names(external.smooths)))
      external.smooths <- external.smooths[c(smth.cols, extra.vars)]
    else
      stop(paste("Do not have the following required columns in data: ", 
                 paste(extra.vars[!(extra.vars %in% names(external.smooths))],collapse=", "), "\n", sep=""))
    
    #Add missing smoothing-parameter columns
    smth.cols <- setdiff(smooth.cols, smth.cols)
    smooth.pars <- rep("Other", length(smth.cols))
    names(smooth.pars) <- smth.cols
    smooth.pars <- rbind(smooth.pars)
    rownames(smooth.pars) <- NULL
    tmp <- cbind(smooth.pars, external.smooths)
    tmp <- tmp[names(smth)]
    #Add extra smooths to smth
    smth <- rbind(smth,tmp)
    #Update the smoothing-parameter schemes
    tmp <- split(tmp, dae::fac.combine(as.list(lapply(tmp[smooth.cols], as.factor)), 
                                       combine.levels = TRUE, sep = "-"))
    sch <- lapply(tmp, function(x) x[1, smooth.cols])
    sch <- do.call(rbind, sch)
    rownames(sch) <- NULL
    spar.schemes <- rbind(spar.schemes, sch[colnames(spar.schemes)])
  } 
  
  #Form a smooths.frame and check that it is valid
  smth <- as.smooths.frame(smth, individuals, times)
  validsmoothsframe <- validSmoothsFrame(smth)  
  if (is.character(validsmoothsframe))
    stop(validsmoothsframe)
  
  #Plot some combination of unsmoothed and smoothed response, AGR and RGR
  if ("profiles" %in% plots)
  { 
    #Plot profiles
    plotSmoothsComparison(data = smth, 
                          response = response, response.smoothed = response.smooth, 
                          times = times, individuals = individuals, 
                          trait.types = traits, 
                          x.title = x.title, y.titles = y.titles, labeller = labeller, 
                          which.plots = "profiles", 
                          plots.include.raw = plots.include.raw,
                          plots.by.pf = plts.by, 
                          facet.x.pf = facet.x.pf, facet.y.pf = facet.y.pf, 
                          colour.pf = colour.pf, colour.column.pf = colour.column.pf, 
                          colour.values.pf = colour.values.pf, 
                          alpha.pf = alpha.pf, 
                          addMediansWhiskers.pf = addMediansWhiskers.pf, 
                          ggplotFuncsProfile = ggplotFuncsProfile, 
                          ggplotFuncsDevnBoxes = ggplotFuncsDevnBoxes, ...)
  }
  
  if ("medians.deviations" %in% plots)
  {
    plotSmoothsMedianDevns(data = smth, 
                           response = response, response.smoothed = response.smooth,
                           times = times, individuals = individuals,
                           x.title = x.title, y.titles = y.titles, labeller = labeller, 
                           plots.by.med = plts.by.med, plots.group.med = plts.group.med,
                           facet.x.med = facet.x.med, 
                           facet.y.med = facet.y.med,
                           trait.types = traits,
                           colour.values.med = colour.values.med, 
                           shape.values.med = shape.values.med, 
                           alpha.med = alpha.med,
                           propn.types.med = propn.types.med, propn.note.med = propn.note.med,
                           ggplotFuncsMedDevn = ggplotFuncsMedDevn, ...)
  }
  
  #Plot some combination of unsmoothed and smoothed response, AGR and RGR
  if (any(c("absolute.boxplots", "relative.boxplots") %in% plots))
  { 
    boxp <- c("absolute.boxplots", "relative.boxplots")[c("absolute.boxplots", "relative.boxplots") %in% plots]
    #Plot boxplots
    plotSmoothsComparison(data = smth, 
                          response = response, response.smoothed = response.smooth, 
                          times = times, individuals = individuals, 
                          trait.types = traits, 
                          plots.by.pf = plts.by, 
                          plots.include.raw = plots.include.raw,
                          x.title = x.title, y.titles = y.titles, labeller = labeller, 
                          which.plots = boxp, 
                          facet.x.pf = facet.x.pf, facet.y.pf = facet.y.pf, 
                          colour.pf = colour.pf, colour.column.pf = colour.column.pf, 
                          colour.values.pf = colour.values.pf, 
                          alpha.pf = alpha.pf, 
                          addMediansWhiskers.pf = addMediansWhiskers.pf, 
                          ggplotFuncsProfile = ggplotFuncsProfile, 
                          ggplotFuncsDevnBoxes = ggplotFuncsDevnBoxes, ...)
  }
  
  invisible(smth)
}

"traitSmooth" <- function(data, response, response.smoothed, individuals, times, 
                          get.rates = TRUE, keep.columns = NULL, 
                          x.title = NULL, y.titles = NULL, 
                          chosen.smooth = list(spline.type = "PS", 
                                               df = NULL, 
                                               lambda = NULL,
                                               smoothing.method = "logarithmic"), 
                          facet.x.chosen = ".", facet.y.chosen = ".", 
                          labeller.chosen = NULL, 
                          colour.chosen = "black", colour.column.chosen = NULL, 
                          colour.values.chosen = NULL, alpha.chosen = 0.3, 
                          addMediansWhiskers.chosen = TRUE, 
                          ggplotFuncsChosen = NULL, 
                          ...)
{
  inargs <- list(...)
  checkEllipsisArgs(c("traitSmooth","probeSmooths","plotProfiles"), inargs)
  
  #Set the list of defaults that are particular to traitSmooth
  probesmth.args <- list(
    smoothing.methods = "logarithmic", 
    spline.types = c("NCSS","PS"), 
    df = 5:7, 
    lambdas = list(PS = round(10^c(-0.5, 0, 0.5, 1), digits = 3)), 
    which.plots = c("profiles", "medians.deviations"), 
    plots.include.raw = TRUE,
    plots.by.pf = "Type", 
    facet.x.pf = c("Method", "Tuning"), 
    facet.x.med = c("Method","Type"), 
    plots.group.med = "Tuning",
    propn.types.med = NULL, 
    printPlot = TRUE)
  
  #Get the arguments that have been set, replace any of the traitSmooth defaults and add extras
  inargs <- list(...)
  if (length(inargs))
    probesmth.args <- utils::modifyList(probesmth.args, inargs, keep.null = TRUE)
  
  tmp <- data
  tmp[times] <- convertTimes2numeric(tmp[[times]])
  #Call probeSmooths
  smth <- do.call(probeSmooths, c(list(data = tmp, 
                                       response = response, response.smoothed, 
                                       times = times, individuals =  individuals, 
                                       get.rates = get.rates, keep.columns = keep.columns,
                                       labeller = labeller.chosen), 
                                  probesmth.args))
  
  #Process chosen model
  if (!is.allnull(chosen.smooth))
  { 
    if (any(unlist(lapply(chosen.smooth, function(x) (length(x) > 1)))))
      stop("the components of chosen.smooth must be single-valued")
    smethods.opt <- c("direct", "logarithmic")
    smethod <- smethods.opt[check.arg.values(chosen.smooth$smoothing.method, options=smethods.opt)]
    #smethods <- smethods.opt[unlist(lapply(probesmth.args$smoothing.methods, check.arg.values, options=smethods.opt))]
    smethods <- levels(smth$Method)
    stypes.opt <- c("NCSS", "PS")
    stype <- stypes.opt[check.arg.values(chosen.smooth$spline.type, options=stypes.opt)]
    #stypes <- stypes.opt[unlist(lapply(probesmth.args$spline.types, check.arg.values, options=stypes.opt))]
    stypes <- levels(smth$Type)
    
    if (!(stype %in% stypes))
      stype <- setdiff(stypes, stype)[1]
    if (!(smethod %in% smethods))
      smethod <- setdiff(smethods, smethod)[1]
    if (is.allnull(chosen.smooth$df) && is.allnull(chosen.smooth$lambda))
    {
      #if (stype == "NCSS"  && !is.allnull(probesmth.args$df))
      #  chosen.smooth$df <- probesmth.args$df[floor((length(probesmth.args$df)+1)/2)]
      if (stype == "NCSS"  && ("df" %in% levels(smth$TunePar)))
      { 
        df <- factor(smth$TuneVal[smth$TunePar == "df"])
        df <- as.numeric(levels(df))
        chosen.smooth$df <- df[floor((length(df)+1)/2)]
      }
      else
      {
        #if (!is.allnull(probesmth.args$lambda))
        #  chosen.smooth$lambda <- 
        #    probesmth.args$lambdas[[stype]][floor((length(probesmth.args$lambdas[[stype]])+1)/2)]
        if ("lambda" %in% levels(smth$TunePar))
        { 
          lambdas <- factor(smth$TuneVal[smth$TunePar == "lambda" & smth$Type == stype])
          lambdas <- as.numeric(levels(lambdas))
          chosen.smooth$lambda <- 
            lambdas[floor((length(lambdas)+1)/2)]
        }
        else
          chosen.smooth$df <- NULL
      }
    }
    if (!is.allnull(chosen.smooth$df))
    {
      tunepar = "df"
      tuneval = as.character(chosen.smooth$df)
    } else
    {
      tunepar = "lambda"
      tuneval = as.character(chosen.smooth$lambda)
    }
    if ((!is.allnull(chosen.smooth$df) && !is.allnull(chosen.smooth$lambda)) ||
        (is.allnull(chosen.smooth$df) && is.allnull(chosen.smooth$lambda)))
      stop("It must be that one and only one of df and lambda is NULL in chosen.smooth")
    combos <- levels(with(smth, fac.combine(list(Type,TunePar,TuneVal,Method), 
                                            combine.levels = TRUE, sep = "-")))
    choice <- paste0(sapply(list(stype, tunepar, tuneval, smethod), as.character), collapse = "-")
    if (!(choice %in% combos))
      stop("The combination of the values of Type, TunePar, TuneVal and Method ", 
           "given in chosen.smooth are not amongst those for the set of smooths that have been investigated")
    
    #Get subset and merge it with the tmp
    smth <- subset(smth, subset = (smth$Type == stype & smth$TunePar == tunepar &
                                     smth$TuneVal == tuneval &  smth$Method == smethod),
                   select = setdiff(names(smth), c("Type","TunePar","TuneVal","Method","Tuning")))
    class(smth) <- "data.frame"
    tmp <- subset(tmp, select = c(individuals, times, setdiff(names(tmp), names(smth))))
    smth <- merge(tmp, smth, sort = FALSE)
    smth <- smth[c(names(data),setdiff(names(smth), names(data)))]
 
    #Plot the profile plots for the chosen smooth
    grates <- NULL
    if (get.rates)
    {
      #Get trait.types set for probeSmooths
      if (!("trait.types" %in% names(probesmth.args)))
        grates <- c("AGR", "RGR")
      else
      {
        options <- c("response", "AGR", "RGR", "all")
        traits <- options[unlist(lapply(probesmth.args$trait.types, check.arg.values, 
                                        options=options))]
        grates <- traits[-match(response, traits)]
        if ("all" %in% traits)
          grates <- c("AGR", "RGR")
        if (!any(c("AGR","RGR") %in% grates))
        {
          warning("trait.types does not include AGR or RGR and so get.rates changed to FALSE")
          grates <- NULL
        }
        else
          grates <- c("AGR","RGR")[c("AGR","RGR") %in% grates]
      }
    }
    if (is.null(grates))
      responses.plot <- response.smoothed
    else
      responses.plot <- c(response.smoothed, paste(response.smoothed, grates, sep = "."))
    
    if (is.null(x.title))
      x.title <- times
    if (is.null(y.titles))
    {
      y.titles <- responses.plot
    }
    else
    {
      if (length(y.titles) != length(responses.plot))
        stop("y.titles is not NULL and have not been provided for the response, the AGR and the RGE")
    }
    names(y.titles) <- responses.plot
    
    #Find out if any plotProfiles arguments that traitSmooths handles have been supplied in '...'
    usedProfile.args <- c("data","times","response","individuals","facet.x","facet.y",
                          "labeller","colour","colour.column","colour.values",
                          "alpha","y.title","printPlot","ggplotFuncs")
    doubleargs <- intersect(names(inargs), usedProfile.args)
    if (length(doubleargs))
      stop("the  'plotProfiles' arguments ",paste0(doubleargs, collapse = ", "),
           " conflict with 'traitSmooth' arguments")
    #extract any remaining plotProfiles arguments from inargs
    pltLongi.args <- setdiff(formalArgs(plotProfiles), c(usedProfile.args))
    pltLongi.args <- names(inargs)[names(inargs) %in% pltLongi.args]
    if (length(pltLongi.args))
      pltLongi.args <- inargs[pltLongi.args]
    else
      pltLongi.args <- NULL
    
    if (is.null(pltLongi.args) || !("title" %in% names(pltLongi.args)))
      title <- paste0("Plot for the choice ", paste0(c(substring(smethod,1,3),stype,
                                                       tunepar,tuneval), collapse = "-"))
    
    #Call plot longitudinal with just the plotLongitudinbal args from inargs
    x.axis <- list(scale_x_continuous(breaks = seq(min(smth[[times]], na.rm = TRUE),
                                                   max(smth[[times]], na.rm = TRUE), 
                                                   by = 2)))
    for (kresp in responses.plot)
      do.call(plotProfiles, c(list(data = smth, times = times,response = kresp,
                                       individuals = individuals,
                                       facet.x = facet.x.chosen,
                                       facet.y = facet.y.chosen,
                                       labeller = labeller.chosen,
                                       colour = colour.chosen,
                                       colour.column = colour.column.chosen,
                                       colour.values = colour.values.chosen,
                                       alpha = alpha.chosen,
                                       title = title,
                                       x.title = x.title,
                                       y.title = y.titles[kresp],
                                       printPlot = TRUE,
                                       ggplotFuncs = c(ggplotFuncsChosen, x.axis)),
                                  pltLongi.args))
  }

  #Make sure that times is of the same type as times in data
  smth[times] <- convertTimesExnumeric(smth[[times]], data[[times]])

  invisible(smth)
}