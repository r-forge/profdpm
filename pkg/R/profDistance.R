profDistance <- function( fit1, fit2 ) {
  ###################################################
  #do some argument checking
  if( !( "profLinear" %in% is(fit1) ) ) { stop("fit1 not of class \'profLinear\'") }
  if( !( "profLinear" %in% is(fit2) ) ) { stop("fit2 not of class \'profLinear\'") }  
  if( length(fit1$clust) != length(fit2$clust) ) { stop("fit1 and fit2 lengths differ") }
  if( length(fit1$group) != length(fit2$group) ) { stop("fit1 and fit2 lengths differ") }
  if( !all.equal(fit1$group, fit2$group) ) { stop("fit1 and fit2 groups differ") }

  uni1 <- unique( cbind( fit1$group, fit1$clust ) )[,2]
  uni2 <- unique( cbind( fit2$group, fit2$clust ) )[,2]
  con <- tot <- 0
  for( i in 1:(length(uni1)-1) ) {
    for( j in (i+1):length(uni2) ) {
      tot <- tot + 1
      if( (uni1[i] == uni1[j] & uni2[i] == uni2[j]) |
          (uni1[i] != uni1[j] & uni2[i] != uni2[j]) ) { con <- con + 1 }
    }
  }
  return( con/tot )
#  ret <- .Call("profDistance", as.integer(uni1[,2]), as.integer(uni2[,2]), PACKAGE="profdpm")
}
