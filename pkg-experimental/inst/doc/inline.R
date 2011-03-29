test <- function() {
    require('inline')
    require('profdpm')
    dll  <- getDynLib("profdpm")
    #regx <- paste(dll[['name']], .Platform$dynlib.ext, '$', sep='')             
    #cat('regex:', regx, '\n')
    #path <- sub(regx, '', dll[['path']])
    #cat('path:', path, '\n')
    #lib  <- paste('-L', path, ' -l', dll[['name']], sep='')
    lib <- paste('-l:', dll[['path']], sep='')
    cat('link:', lib, '\n')
    sig  <- signature()
    cat('sig:', sig, '\n')
    bod  <- 'printf("Hello World!\\n");\nreturn R_NilValue;'
    test <- cfunction(sig, bod, libargs=lib, language='C', verbose=TRUE)
    return(test)
}
    
