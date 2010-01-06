profLinear <- function(y, x, group, clust, param, method="Shotwell", stop=FALSE,
                       maxiter=1000, crit=1e-5, prior="Dirichlet", verbose=FALSE) {
  ###################################################
  #do some argument checking
  if(!is.numeric(y)) { stop("y must be numeric") }
  if(!is.matrix(x)|!is.numeric(x)) { stop("x must be a numeric matrix") }
  if(missing(group)) { group <- seq(1, length(y)) }
  if(length(y) != length(group)) { stop("length(y) must equal length(group)") }
  if(length(y) != nrow(x)) { stop("length(y) must equal nrow(x)") }
  if(missing(param)) { param <- list(alpha=1,a0=0.001,b0=0.001,m0=rep(0,ncol(x)),s0=1.000) }
  if(missing(clust)) { clust <- FALSE }
  else { if(length(y) != length(clust)) { stop("length(y) must equal length(clust)") } }
  if(!is.numeric(stop) & !(is.logical(stop) & stop==FALSE) ) { stop("stop must be FALSE or numeric") }

  ###################################################
  #remove missing observations, issue warning
  miss <- apply( is.na( cbind( y, x ) ), 1, any ) 
  ry <- y[!miss]
  rx <- x[!miss,]
  rg <- group[!miss]
  if( is.logical(clust) ) { rc <- FALSE }
  else{ rc <- clust[!miss] }
  if( any( miss ) ) {
    warning( "removed observations with missing values: ", paste(" ", (1:length(y))[miss], sep="") )
  }

  ###################################################
  #order the data according to group
  #convert ordered y to double
  #convert ordered group to integers from 0,1,...
  #convert ordered clust to integers from 0,1,...
  rg <- factor(rg)
  ord <- order(rg)
  ry <- as.double(ry[ord])
  rg <- as.integer(unclass(rg[ord])-1)
  if( !is.logical(rc) ) { rc <- as.integer(unclass(rc[ord])-1) }

  ###################################################
  #transpose ordered x to simplify BLAS/LAPACK calls
  #matrices are double storage by column
  rx <- t(as.matrix(x[ord,]))

  ###################################################
  #convert prior to integer
  if( prior == "Dirichlet" ) { prior <- 0 }
  else if( prior == "cluster" ) { prior <- 1 }
  else{ 
    prior <- 0
    warning("prior must be \'Dirichlet\' or \'cluster\', using \'Dirichlet\'")  
  }

  ###################################################
  #convert method to integer
  if( method == "none" ) { method <- 0 }
  else if( method == "Shotwell" )      { method <- 1 }
  else if( method == "agglomerative" ) { method <- 2 }
  else {
    method <- 1 #default is "Shotwell"
    warning("method must be \'Shotwell\', \'agglomerative\', or \'none\'", )
  }

  ###################################################
  #convert stop to integer
  if( is.logical(stop) ) { stop <- -1 }

  ###################################################
  #call the C function
  ret <- .Call("profLinear", ry, rx, rg, rc, as.list(param), as.integer(method), as.integer(stop), 
                as.integer(maxiter), as.numeric(crit), as.integer(prior), as.logical(verbose), PACKAGE="profdpm")

  ###################################################
  #undo ordering
  ret$y[ord] <- ret$y
  ret$x <- t(ret$x)
  ret$x[ord,] <- ret$x
  ret$group[ord] <- ret$group
  ret$clust[ord] <- ret$clust
  ret$clust <- unclass(as.factor(ret$clust))
  attributes(ret$clust) <- NULL
  return(ret)  
}


