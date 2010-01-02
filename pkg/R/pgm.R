pgm <- function( filename ) {
  ###################################################
  ret <- .Call("pgm", as.character(filename), PACKAGE="profdpm")
  return( ret )
}
