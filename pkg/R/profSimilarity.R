profSimilarity <- function( obj1, obj2 ) {
  ###################################################
  #do some argument checking
  if( !( "profLinear" %in% is(obj1) ) ) { stop("obj1 not of class \'profLinear\'") }
  if( !( "profLinear" %in% is(obj2) ) ) { stop("obj2 not of class \'profLinear\'") }  
  if( length(obj1$clust) != length(obj2$clust) ) { stop("obj1 and obj2 lengths differ") }
  if( length(obj1$group) != length(obj2$group) ) { stop("obj1 and obj2 lengths differ") }
  if( !all.equal(obj1$group, obj2$group) ) { stop("obj1 and obj2 groups differ") }

  uni1 <- unique( cbind( obj1$group, obj1$clust ) )[,2]
  uni2 <- unique( cbind( obj2$group, obj2$clust ) )[,2]
  ret <- .Call("profSimilarity", as.integer(uni1), as.integer(uni2), PACKAGE="profdpm")
  return( ret )
}
