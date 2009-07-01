profLinear <- function(form, data, group, parm, iter=1000, crit=0.001) {
  if(missing(group)) { group <- 0:nrow(data) }
  .Call("profLinear", as.formula(form), as.data.frame(data), as.integer(group), as.list(parm), as.integer(iter), as.numeric(crit), PACKAGE="profdpm")
} 
