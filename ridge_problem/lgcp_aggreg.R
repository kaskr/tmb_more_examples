library(TMB)
compile("lgcp_aggreg.cpp","-O0 -g")

## nobs <- 100
## ngroup <- 3

nobs <- 1
ngroup <- 2

map <- list(mu=factor(c(1,NA)))
parameters <- list(
                   mu = c(-2,-2),
                   logsd = 0,
                   x=array(0,c(ngroup,nobs))
                   )
set.seed(123)
## N <- replicate(nobs,
##                rpois(1,sum(exp(rnorm(ngroup,parameters$mu,exp(parameters$logsd)))))
##                )
N <- 2
mu <- -2
Q <- 1
N <- 2*(exp(mu+1)+Q)

data <- list(N=N)

dyn.load("lgcp_aggreg.so")
obj <- MakeADFun(data,parameters,random="x",map=map)
obj$env$random.start <- expression(c(par.fixed[1],1))


####################################################################
newmyf <- function(mu){
    grid <- seq(-3,4.5,by=0.05) - 2
    df <- expand.grid(grid,grid)
    df <- cbind(mu,0,df)
    mat <- t(as.matrix(df))
    val <- apply(mat,2,obj$env$f)
    dim(val) <- rep(length(grid),2)
    image(grid,grid,exp(-val))
    contour(grid,grid,exp(-val),add=TRUE)
    f <- function(eta)log(N-exp(eta))
    plot(f,-5,2,add=TRUE,n=1000,lwd=2)
    lap <- exp(-obj$fn(c(mu,0)))
    p <- obj$env$last.par[-(1:2)]
    abline(0,1,col="green",lty="dashed")
    points(p[1],p[2],pch=16)
    prior <- obj$env$parList()$mu
    points(prior[1],prior[2],pch=16,col="blue")
    euler <- sum(exp(-val))*(grid[2]-grid[1])^2
    txt <- paste("euler:",round(euler,6),"laplace:",round(lap,6),"mu1:",mu)
    title(txt)
    print(txt)
}
pdf("surface30.pdf")
gr <- seq(-.1,.1,by=0.01) + -2
for(mu in gr){
    newmyf(mu)
}
dev.off()


