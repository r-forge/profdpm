summary.profBinary <-
function (x, ...) 
{
    r <- length(unique(x$clust))
    q <- ncol(x$y)
    tclust <- table(x$clust)
    term.labels <- attr(terms(x$model), "term.labels")
    if (attr(terms(x$model), "intercept")) 
        term.labels <- c("(intercept)", term.labels)
    cat("----------\n")
    out <- vector("list", length = r)
    for (cls in 1:r) {
        out[[cls]]$groups <- tclust[cls]
        cat("cluster: ", cls, "\n", sep = "")
        cat("groups:  ", tclust[cls], "\n", sep = "")
        M <- x$a[[cls]] / (x$a[[cls]] + x$b[[cls]])
        lower <- qbeta(0.025, x$a[[cls]], x$b[[cls]])
        upper <- qbeta(0.975, x$a[[cls]], x$b[[cls]])
        info <- data.frame(estimate = M, lower95 = lower, 
            upper95 = upper)
        row.names(info) <- term.labels
        out[[cls]]$summary <- info
        print(format(info, nsmall = 3, digits = 3))
        cat("----------\n")
    }
    invisible(out)
}

