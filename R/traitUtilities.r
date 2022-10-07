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

#Function to combine just the smoothing-parameter factors 
fac.tunecombine <- function(dat, factors, smooth.cols)
{
  tune.facs <- intersect(factors, smooth.cols)
  if ("Method" %in% tune.facs)
    dat["Method"] <- dae::fac.recast(dat[["Method"]], 
                                     newlevels = substring(levels(dat[["Method"]]), 
                                                           first = 1, last = 3))
  if (length(tune.facs) <= 0) 
    fac <- NULL  
  else
  {
    if (length(tune.facs) == 1)
      fac <- dat[[tune.facs]]
    else
      fac <- with(dat, dae::fac.combine(as.list(dat[tune.facs]),
                                        combine.levels = TRUE, sep = "-"))
  }
  return(fac)
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

#Function to set up a facet using fac.mixcombine, fac.tunecombine and fac.combine
"setupFacet" <- function(data, facet, collapse.facets, combined.name = "SmoothParams", smooth.cols)
{
  newfacet <- facet
  if (any(smooth.cols %in% newfacet))
  {
    if (collapse.facets)
    {  
      data[combined.name] <- fac.mixcombine(data, facet, smooth.cols = smooth.cols)
      newfacet <- combined.name
    } else
    {
      data[combined.name] <- fac.tunecombine(data, facet, smooth.cols = smooth.cols)
      newfacet <- c(setdiff(newfacet, smooth.cols), combined.name)
    }
  } else #no smooth.cols & several factors to collapse
  { 
    if (!all(newfacet == ".") && (length(newfacet) > 1) && collapse.facets)
    { 
      data[combined.name] <- dae::fac.combine(as.list(data[newfacet]), combine.levels = TRUE)
      newfacet <- combined.name
    }
  }
  return(list(newfacet = newfacet, data = data))
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
    stop(paste0("\n The number of different combinations of (i) the smoothing-parameter values that ", 
                "are available and (ii) the levels combination of the following factors nominated ",
                "in the facet/plots arguments are not equal: ", paste0(spar.facs, collapse = ", ")))
  return(data)
}  

#Function to set up scale_x_continuous for times
setScaleTime <- function(times, breaks.spacing.x = -2)
{
  abs.spacing.x <- abs(breaks.spacing.x)
  time.vals <- sort(unique(times))
  time.vals <- time.vals[!is.na(time.vals)]
  brks <- seq(min(times, na.rm = TRUE),
              max(times, na.rm = TRUE), 
              by = abs.spacing.x)
  if (breaks.spacing.x < 0)
    brks <- brks[brks %in% time.vals]
  minbrks <- seq(min(time.vals), max(time.vals), 
                 by = abs.spacing.x/2)
  if (breaks.spacing.x < 0)
  {  
    minbrks <- minbrks[minbrks %in% time.vals]
    if (length(minbrks) == length(brks) && all(minbrks == brks))
    { 
      minbrks <- seq(min(time.vals), max(time.vals), 
                     by = 1)
      if (breaks.spacing.x < 0)
        minbrks <- minbrks[minbrks %in% time.vals]
    }
  }
  scale.time <- scale_x_continuous(limits = range(time.vals),
                                   breaks = brks, minor_breaks = minbrks)
  return(scale.time)
}

#Function to put options for median deviations plots into a list
args4meddevn.plot <- function(plots.by = NULL, plots.group = "Tuning", 
                              facet.x = c("Method","Type"), facet.y = ".", 
                              facet.labeller = NULL, facet.scales = "free_x",
                              breaks.spacing.x = -4, angle.x = 0, 
                              colour.values = NULL, shape.values = NULL, alpha = 0.5, 
                              propn.note = TRUE, propn.types = NULL, 
                              ggplotFuncs = NULL, 
                              ...)
{
  inargs <- list(...)
  checkEllipsisArgs("args4meddevn.plot", inargs)
  
  #Checking of the arguments that control the plots layout for the medians deviations plots 
  checkLayoutArgs(data = NULL, plots.by, plts.group = plots.group, facet.x, facet.y)
  
  med.opts <- list(plots.by = plots.by, 
                   plots.group = plots.group,
                   facet.x = facet.x, facet.y = facet.y, 
                   facet.labeller = facet.labeller, 
                   facet.scales = facet.scales, 
                   breaks.spacing.x = breaks.spacing.x,
                   angle.x = angle.x, 
                   colour.values = colour.values, 
                   shape.values = shape.values, alpha = alpha, 
                   propn.note = propn.note, propn.types = propn.types, 
                   ggplotFuncs = ggplotFuncs)
  return(med.opts)
}

#Function to put options for the chosen plots into a list
args4chosen.plot <- function(plots.by = NULL, 
                              facet.x = ".", facet.y = ".", 
                              include.raw = "no", 
                              collapse.facets.x = FALSE, collapse.facets.y = FALSE, 
                              facet.labeller = NULL, facet.scales = "fixed", 
                              breaks.spacing.x = -2, angle.x = 0, 
                              colour = "black", colour.column = NULL, 
                              colour.values = NULL, alpha = 0.3, 
                              addMediansWhiskers = TRUE,
                              ggplotFuncs = NULL, 
                              ...)
{
  inargs <- list(...)
  checkEllipsisArgs("args4chosen.plot", inargs)
  
  pf.opts <- list(plots.by = plots.by, 
                  facet.x = facet.x, facet.y = facet.y, 
                  include.raw = include.raw, 
                  collapse.facets.x = collapse.facets.x, 
                  collapse.facets.y = collapse.facets.y, 
                  facet.labeller = facet.labeller, 
                  facet.scales = facet.scales, 
                  breaks.spacing.x = breaks.spacing.x, 
                  angle.x = angle.x, 
                  colour = colour, colour.column = colour.column, 
                  colour.values = colour.values, alpha = alpha, 
                  addMediansWhiskers = addMediansWhiskers,
                  ggplotFuncs = ggplotFuncs)
  return(pf.opts)
}

#Function to put smoothing parameters options into a list
args4chosen.smooth <- function(smoothing.methods = "logarithmic", 
                               spline.types = "PS", 
                               df = NULL, 
                               lambdas = NULL,
                               combinations = "single")
{
  options <- c("single")
  comb.opt <- options[check.arg.values(combinations, options=options)]
  if (comb.opt != "single")
    stop("combinations must be set to single for chosen.smooth")
  
  #Deal with smoothing arguments
  smethods.opt <- c("direct", "logarithmic")
  smethods <- smethods.opt[unlist(lapply(smoothing.methods, check.arg.values, options=smethods.opt))]
  stypes.opt <- c("NCSS", "PS")
  stypes <- stypes.opt[unlist(lapply(spline.types, check.arg.values, options=stypes.opt))]
  options <- c("differences","derivatives")
  
  smth.params <- list(smoothing.methods = smethods, 
                      spline.types = stypes, 
                      df = df, 
                      lambdas = lambdas, 
                      combinations = combinations)
  return(smth.params)
}

#Function to put options for profile plots into a list
args4profile.plot <- function(plots.by = "Type", 
                              facet.x = c("Method", "Tuning"), facet.y = ".", 
                              include.raw = "facet.x", 
                              collapse.facets.x = TRUE, collapse.facets.y = FALSE, 
                              facet.labeller = NULL, facet.scales = "fixed", 
                              breaks.spacing.x = -4, angle.x = 0, 
                              colour = "black", colour.column = NULL, 
                              colour.values = NULL, alpha = 0.3, 
                              addMediansWhiskers = TRUE,
                              ggplotFuncs = NULL, 
                              ...)
{
  inargs <- list(...)
  checkEllipsisArgs("args4profile.plot", inargs)
  
  checkLayoutArgs(data = NULL, plots.by, plts.group = NULL, facet.x, facet.y)
  
  pf.opts <- list(plots.by = plots.by, 
                  facet.x = facet.x, facet.y = facet.y, 
                  include.raw = include.raw, 
                  collapse.facets.x = collapse.facets.x, 
                  collapse.facets.y = collapse.facets.y, 
                  facet.labeller = facet.labeller, 
                  facet.scales = facet.scales, 
                  breaks.spacing.x = breaks.spacing.x, 
                  angle.x = angle.x, 
                  colour = colour, colour.column = colour.column, 
                  colour.values = colour.values, alpha = alpha, 
                  addMediansWhiskers = addMediansWhiskers,
                  ggplotFuncs = ggplotFuncs)
  return(pf.opts)
}

#Function to put smoothing parameters options into a list
args4smoothing <- function(smoothing.methods = "logarithmic", 
                           spline.types = c("NCSS","PS"), 
                           df = 5:7, 
                           lambdas = list(PS = round(10^c(-0.5, 0, 0.5, 1), 
                                                     digits = 3)),
                           smoothing.segments = NULL, npspline.segments = NULL, 
                           na.x.action="exclude", na.y.action = "trimx", 
                           external.smooths = NULL, 
                           correctBoundaries = FALSE, 
                           combinations = "allvalid", ...)
{
  inargs <- list(...)
  checkEllipsisArgs("args4smoothing", inargs)
  
  options <- c("allvalid", "parallel", "single")
  #Deal with smoothing arguments
  smethods.opt <- c("direct", "logarithmic")
  smethods <- smethods.opt[unlist(lapply(smoothing.methods, check.arg.values, options=smethods.opt))]
  stypes.opt <- c("NCSS", "PS")
  stypes <- stypes.opt[unlist(lapply(spline.types, check.arg.values, options=stypes.opt))]
  options <- c("differences","derivatives")
  
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
  
  smth.params <- list(smoothing.methods = smethods, 
                      spline.types = stypes, df = df, lambdas = lambdas, 
                      smoothing.segments = smoothing.segments, 
                      npspline.segments = npspline.segments, 
                      na.x.action = na.x.action, na.y.action = na.y.action, 
                      external.smooths = external.smooths, 
                      correctBoundaries = correctBoundaries, 
                      combinations = combinations)
  return(smth.params)
}

conv2TuneParams <- function(smth.args, smth)
{
  smethods.opt <- c("direct", "logarithmic")
  smethod <- smethods.opt[check.arg.values(smth.args$smoothing.method, options=smethods.opt)]
  smethods <- levels(smth$Method)
  stypes.opt <- c("NCSS", "PS")
  stype <- stypes.opt[check.arg.values(smth.args$spline.type, options=stypes.opt)]
  stypes <- levels(smth$Type)
  
  if (!(stype %in% stypes))
    stype <- setdiff(stypes, stype)[1]
  if (!(smethod %in% smethods))
    smethod <- setdiff(smethods, smethod)[1]
  if (is.allnull(smth.args$df) && is.allnull(smth.args$lambdas))
  {
    if (stype == "NCSS"  && ("df" %in% levels(smth$TunePar)))
    { 
      df <- factor(smth$TuneVal[smth$TunePar == "df"])
      df <- as.numeric(levels(df))
      smth.args$df <- df[floor((length(df)+1)/2)]
    }
    else
    {
      if ("lambda" %in% levels(smth$TunePar))
      { 
        lambdas <- factor(smth$TuneVal[smth$TunePar == "lambda" & smth$Type == stype])
        lambdas <- as.numeric(levels(lambdas))
        smth.args$lambdas <- 
          lambdas[floor((length(lambdas)+1)/2)]
      }
      else
        smth.args$df <- NULL
    }
  }
  if (!is.allnull(smth.args$df))
  {
    tunepar = "df"
    tuneval = as.character(smth.args$df)
  } else
  {
    tunepar = "lambda"
    tuneval = as.character(smth.args$lambdas)
  }
  return(list(stype = stype, tunepar = tunepar, tuneval = tuneval, smethod = smethod))
}

#Function to generate the set of smoothing schemes
makeSmoothSchemes <- function(combinations, smethods, stypes, df, lambdas)
{
  #Construct combinations of smoothing parameters
  if (combinations == "allvalid")
  {
    if (!is.allnull(lambdas))
    {  
      if  (!is.list(lambdas))
      {  
        #warning("lambdas have been converted to a list with the single component names PS")
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
          stop("one of df or lambdas with NCSS must be specified for spline.type NCSS")
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
            stop("At least one value needs to be specfied for lambdas when spline.type is PS")
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
    if (combinations == "parallel")
    { 
      #Check all f the same length
      if (!all(sapply(list(smethods, stypes, df, lambdas), 
                      function(x, first) length(x) == first, first = length(smethods))))
        stop(paste0("For combinations set to parallel, smoothing mewthods, spline.types, df and lambdas must be of the same length ", 
                    "(if df or lambdas not NULL pad with NAs)"))
      
      TVal <- df
      TVal[which(!is.na(lambdas))] <- lambdas[which(!is.na(lambdas))]
      TPar <- rep("df", length(df))
      TPar[which(!is.na(lambdas))] <- "lambda"
      
      spar.schemes <- data.frame(Type = stypes,
                                 TunePar = TPar, 
                                 TuneVal = TVal, 
                                 Tuning = paste(TPar, TVal, sep = "-"), 
                                 Method = smethods)
      attr(spar.schemes, which = "nschemes") <- nrow(spar.schemes)
    } else
      stop("The setting of combinations in smoothing parameters must be either allvalid or parallel")
  }
  return(spar.schemes) 
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
    if (length(grates) > 0)
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

