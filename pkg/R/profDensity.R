profDensity <- function(data, parm=list(alpha=0.001,a0=0.001,b0=0.001,m0=0,s0=1),iter=1000,crit=0.001) {
  .Call("profDensity", as.numeric(data), as.list(parm),as.integer(iter), as.numeric(crit), PACKAGE="profdpm")
}

