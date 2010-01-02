Rand <- function( x1, x2 ) {
  ###################################################
  #do some argument checking
  if( length( x1) != length( x2 ) ) { stop("x1 and x2 lengths differ") }
  ret <- .Call("Rand", as.integer(x2), as.integer(x1), PACKAGE="profdpm")
  return( ret )
}
