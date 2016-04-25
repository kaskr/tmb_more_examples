## Simulate data
source("simdata.R")

## Adapt parameter list
parameters$u <- rep(0,n)
parameters$uaux <- rep(0,n)

require(TMB)
compile('transform.cpp')
dyn.load(dynlib('transform'))

data$do_transform <- 1
data$tiny <- 1e-2
obj <- MakeADFun(data, parameters, random=c("u","uaux"), DLL="transform")

obj$fn()

system.time( opt <- nlminb(obj$par,obj$fn,obj$gr) )
(sdr <- sdreport(obj))
