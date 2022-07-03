#A function to compute a FUN from the response for each individual
"byIndv_ValueCalc" <- function(data, response, individuals = "Snapshot.ID.Tag", 
                               FUN = "max", which.obs = FALSE, which.values = NULL, 
                               addFUN2name = TRUE, 
                               weights=NULL, na.rm=TRUE, sep=".", ...)
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
  #Check that response, individuals, which.values and weights are in data
  Names <- c(response, individuals, which.values, weights)
  checkNamesInData(Names, data = data)

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
                        if (all(is.na(data[[response]])))
                          vals <- NA
                        else
                          vals <- FUNC(x=data[[response]], na.rm=na.rm, ...) 
                        return(vals)
                      },
                      response=response, FUNC=funct, na.rm=na.rm, ...)
  else
    val.dat <- lapply(data, 
                      function(data, response, weights, FUNC, na.rm, ...)
                      { 
                        if (all(is.na(data[[response]])))
                          vals <- NA
                        else
                          vals <- FUNC(x=data[[response]], w = data[[weights]], na.rm=na.rm, ...) 
                        return(vals)
                      },
                      response=response, weights=weights, FUNC=funct, na.rm=na.rm, ...)
  val.dat <- as.data.frame(do.call(rbind, val.dat))
  if (addFUN2name)
    names(val.dat) <- paste(response,FUN,sep=".")
  else
    names(val.dat) <- response
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
  if (which.obs || !is.null(which.values))
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
                            if (is.na(w))
                              val <- NA
                            else
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
      {
        resp.which <- paste(resp.val,which.values,sep=".") 
        which.dat[1] <- convertTimes2numeric(which.dat[[1]])
      }
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
