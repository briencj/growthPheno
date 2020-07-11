"twoLevelOpcreate" <- function(responses, data, treatment.factor="Treatment.1", 
                        suffices.treatment=c("Cont","Salt"), control=1, 
                        columns.suffixed = NULL, 
                        operations = "/", suffices.results="OST", 
                        columns.retained = c("Snapshot.ID.Tag",
                                             "Smarthouse","Lane", 
                                             "Zones","xZones","SHZones",
                                             "ZLane","ZMainplots","xMainPosn", 
                                             "Genotype.ID"),
                        by = c("Smarthouse","Zones","ZMainplots"))
{ 
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