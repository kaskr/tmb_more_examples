library(TMB)

## compile("hmm.cpp",CPPFLAGS="-DEIGEN_USE_BLAS")
compile("hmm.cpp")
dyn.load(dynlib("hmm"))

## Simulate data
source("simdata.R")

## Adapt data for HMM
data$grid <- seq(0,1,length=101)

obj <- MakeADFun(data=data,parameters=parameters,DLL="hmm")
obj$fn()

system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
(sdr <- sdreport(obj,opt$par))
