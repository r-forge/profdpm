profLinear <- function(y, x, group, parm, iter=1000, crit=0.001) {
  ###################################################
  #do some argument checking
  if(!is.numeric(y)) { stop("y must be numeric") }
  if(!is.matrix(x)|!is.numeric(x)) { stop("x must be a numeric matrix") }
  if(missing(group)) { group <- seq(1, length(y)) }
  if(length(y) != length(group)) { stop("length(y) must equal length(group)") }
  if(length(y) != nrow(x)) { stop("length(y) must equal nrow(x)") }

  ###################################################
  #order the data according to group
  #convert ordered y to double
  group <- as.factor(group)
  ord <- order(group)
  ry <- as.double(y[ord])

  ###################################################
  #transpose ordered x to simplify BLAS/LAPACK calls
  #matrices are double storage by column
  rx <- t(as.matrix(x[ord,]))

  ###################################################
  #convert ordered group to integers from 0,1,...
  rg <- as.integer(unclass(group[ord])-1)
 
  ###################################################
  #call the C function
  ret <- .Call("profLinear", ry, rx, rg, as.list(parm), as.integer(iter), as.numeric(crit), PACKAGE="profdpm")

  ###################################################
  #undo ordering
  ret$y[ord] <- ret$y
  ret$x <- t(ret$x)
  ret$x[ord,] <- ret$x
  ret$group[ord] <- ret$group
  ret$clust[ord] <- ret$clust
  return(ret)  
}


