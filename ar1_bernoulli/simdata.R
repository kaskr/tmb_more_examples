## Simulate data
set.seed(123)
n <- 10000
sigma <- .3
phi <- .8
simdata <- function(){
    u <- numeric(n)
    u[1] = rnorm(1)
    if(n>=2)
        for(i in 2:n){
            u[i] = phi * u[i-1] + rnorm(1, sd = sqrt(1 - phi^2))
        }
    u <- u * sigma
    x <- as.numeric( rbinom(n, 1, plogis(u)) )
    data <- list(obs=x)
    data
}
##
data <- simdata()
parameters <- list(phi=phi, logSigma=log(sigma))
