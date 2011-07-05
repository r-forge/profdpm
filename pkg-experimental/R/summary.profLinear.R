summary.profLinear <-
function(x, ...) {
    r <- length(unique(x$clust))
    q <- length(x$m[[1]])
    uclust <- x$clust[match(unique(x$group), x$group)]
    tclust <- table(uclust)
    tclusto <- table(x$clust)
    term.labels <- attr(terms(x$model), 'term.labels')
    if(attr(terms(x$model), 'intercept'))
        term.labels <- c('(intercept)', term.labels)
    cat('----------\n')
    out <- vector('list', length=r)
    for(cls in 1:r) {
        out[[cls]]$groups <- tclust[cls]
        out[[cls]]$obs <- tclusto[cls]
        cat('cluster: ', cls, '\n', sep='')        
        cat('groups:  ', tclust[cls], '\n', sep='') 
        cat('obs:     ', tclusto[cls], '\n', sep='')
        #Laplace approximation to 95%CI for coefficients
        M <- x$m[[cls]]
        S <- 1/sqrt(diag(x$s[[cls]] * (x$a[[cls]]+q/2) / x$b[[cls]])) 
        lower <- M - 1.92 * S
        upper <- M + 1.92 * S
        info <- data.frame(estimate=x$m[[cls]], 
            lower95=lower, upper95=upper)
        row.names(info) <- term.labels
        out[[cls]]$summary <- info
        print(format(info, nsmall=3, digits=3))
        cat('----------\n')
    }
    invisible(out)
}

