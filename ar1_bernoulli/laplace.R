## Simulate data
source("simdata.R")

## Adapt parameter list
parameters$u <- rep(0,n)

require(TMB)
compile('laplace.cpp')
dyn.load(dynlib('laplace'))

obj <- MakeADFun(data, parameters, random="u", DLL="laplace")
obj$fn()

system.time( opt <- nlminb(obj$par,obj$fn,obj$gr) )
(sdr <- sdreport(obj))
