plot.profDensity <- function(x, index, ...) {
  if(missing(index)) {
    hist(x$data, freq=FALSE, main="", xlab="data", ...)
    support <- seq(min(x$data), max(x$data), length.out=1000)
    density <- rep(0, 1000)
    counts <- table(x$index)
    for(i in 1:length(counts)) {
      p <- counts[i] / sum(counts)
      mu <- x$m[i]
      sd <- sqrt(x$b[i] / x$a[i])
      density <- density + p*dnorm(support,mu,sd)
    }
  }
  else {
    sub <- subset(x$data, subset = (x$index == index))
    hist(sub, freq = FALSE, main="", xlab=paste("data group ",index,sep=""), ...)
    support <- seq(min(sub), max(sub), length.out=1000)
    density <- dnorm(support, x$m[index], sqrt(x$b[index]/x$a[index]))
  }
  lines(support, density, col="red")  
}

