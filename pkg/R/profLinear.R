profLinear <- function(y, x, group, parm, iter=1000, crit=1e-5, verb=FALSE) {
  ###################################################
  #do some argument checking
  if(!is.numeric(y)) { stop("y must be numeric") }
  if(!is.matrix(x)|!is.numeric(x)) { stop("x must be a numeric matrix") }
  if(missing(group)) { group <- seq(1, length(y)) }
  if(length(y) != length(group)) { stop("length(y) must equal length(group)") }
  if(length(y) != nrow(x)) { stop("length(y) must equal nrow(x)") }
  if(missing(parm)) { parm <- list(alpha=0.001,a0=0.001,b0=0.001,m0=rep(0,ncol(x)),s0=1.000) }

  ###################################################
  #remove missing observations, issue warning
  miss <- apply( is.na( cbind( y, x ) ), 1, any ) 
  ry <- y[!miss]
  rx <- x[!miss,]
  rg <- group[!miss]
  if( any( miss ) ) {
    warning( "removed observations with missing values: ", (1:length(y))[miss] )
  }

  ###################################################
  #order the data according to group
  #convert ordered y to double
  rg <- as.factor(rg)
  ord <- order(rg)
  ry <- as.double(ry[ord])

  ###################################################
  #transpose ordered x to simplify BLAS/LAPACK calls
  #matrices are double storage by column
  rx <- t(as.matrix(x[ord,]))

  ###################################################
  #convert ordered group to integers from 0,1,...
  rg <- as.integer(unclass(rg[ord])-1)
 
  ###################################################
  #call the C function
  ret <- .Call("profLinear", ry, rx, rg, as.list(parm), as.integer(iter), as.numeric(crit), as.logical(verb), PACKAGE="profdpm")

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


