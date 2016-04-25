## Simulate data
source("simdata.R")

## Adapt parameter list
parameters$u <- rep(0,n)
parameters$utilde <- rep(0,n)

require(TMB)
compile('test.cpp')
dyn.load(dynlib('test'))

data$do_transform <- 1
data$tiny <- 1e-2
obj <- MakeADFun(data, parameters, profile=c("u"), DLL="test")

obj$fn()

system.time( opt <- nlminb(obj$par,obj$fn,obj$gr) )
(sdr <- sdreport(obj))
