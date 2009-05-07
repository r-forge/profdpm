profLinear <- function(form, data, parm, iter=1000, crit=0.001) {
  .Call("profLinear", as.formula(form), as.data.frame(data), as.list(parm), as.integer(iter), as.numeric(crit), PACKAGE="profdpm")
} 
