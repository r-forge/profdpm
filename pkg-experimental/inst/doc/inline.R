profBinaryExtension <- function() {
    require('inline')
    require('profdpm')
    dll  <- getDynLib("profdpm")
    inc  <- 
    lib  <- paste('-l:', dll[['path']], sep='')
    sig  <- signature(x='numeric')
    bod  <- "
        if(!isNumeric(x) || length(x) < 1)
            error(\"invalid %s argument\", \"'x'\");
        printf(\"first value: %f\\n\", REAL(x)[0]);
        return R_NilValue;
    "
    profBinary  <- cfunction(sig, bod, libargs=lib, verbose=TRUE)
    function(x) invisible(profBinary(x))
}    

f <- profBinaryExtension()
