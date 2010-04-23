profBinary <- function(y, clust, param, method="stochastic", 
                       maxiter=1000, crit=1e-5, verbose=FALSE) {
  ###################################################
  #do some argument checking
  if(!is.numeric(y)) { 
    stop("y must be numeric") 
  }
  y <- as.matrix(y)
  if(missing(clust)) { 
    clust <- FALSE 
  }
  else { 
    clust <- as.factor(clust)
    if(nrow(y) != length(clust)) { 
      stop("nrow(y) must equal length(clust)") 
    } 
  }
  if(missing(param)) { 
    param <- list(alpha=1,a0=1.00,b0=1.00) 
  }
  else if(!is.list(param)) {
    warning("param must be a list, using defaults")
    param <- list(alpha=1,a0=1.00,b0=1.00) 
  } 
  else{
    if(length(names(param)) == 0) {
      warning("param argument does not include any named items, using defaults")
      param <- list(alpha=1,a0=1.00,b0=1.00)
    }
    else if(length(param) > length(names(param))) {
      warning("param contains unnamed items")
    }
  }
  if(!is.character(method)) {
    warning("method must be a character string, using default")
    method <- "stochastic"
  }
  if(!is.numeric(maxiter) | maxiter < 0) {
    warning("maxiter must be numeric and non-negative, using default")
    maxiter <- 1000
  }
  if(!is.numeric(crit) | crit < 0) { 
    warning("crit must be numeric and non-negative, using default") 
    crit <- 1e-5
  }
  if(!is.logical(verbose)) {
    warning("verbose must be a logical, using default")
    verbose <- FALSE
  }

  ###################################################
  #remove missing observations, issue warning
  miss <- apply( is.na( y ), 1, any ) 
  ry <- y[!miss,]
  if( is.logical(clust) ) { rc <- FALSE }
  else{ rc <- unclass(clust[!miss])-1 }
  if( any( miss ) ) {
    warning( "removed observations with missing values: ", 
      paste(" ", which(miss), sep="") )
  }

  ###################################################
  #convert method to integer
  if(      method == "none" )          { method <- 0 }
  else if( method == "stochastic" )    { method <- 1 }
  else if( method == "agglomerative" ) { method <- 2 }
  else if( method == "gibbs" )         { method <- 3 }
  else {
    method <- 1 #default is "stochastic"
    warning("method must be \'stochastic\', \'agglomerative\', or \'none\'", )
  }

  ###################################################
  #convert ry to integer storage
  storage.mode(ry) <- "integer"

  ###################################################
  #call the C function
  ret <- .Call("profBinary", ry, rc, as.list(param), as.integer(method),
                as.integer(maxiter), as.double(crit), as.logical(verbose), 
                PACKAGE="profdpm")

  ###################################################
  #adjust clust
  #ret$clust <- unclass(as.factor(ret$clust))
  #attributes(ret$clust) <- NULL
  return(ret)  
}
