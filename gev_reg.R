rm(list=ls())
gc()
library(evd)
## function 
gevreg = function(x, z, lambda = 0)
{
  
  l2gev = function (tvec, lambda)
  {
    loc.vec = zz%*%tvec[1:(p+1)]
    # loglikelihood
    v1 = - sum(lgev(x, loc = loc.vec, scale = tvec[p+2], shape = tvec[p+3]))
    # regularization
    v2 = lambda*sum(tvec[2:(p+1)]^2)
    v = v1 + v2
    return(v)
  }
  
  lgev = function (x, loc = 0, scale = 1, shape = 0) 
  {
    if (min(scale) <= 0) 
      return( - 1e+6)
    if (length(shape) != 1) 
      stop("invalid shape")
    x <- (x - loc)/scale
    if (shape == 0) 
      d <- log(1/scale) - x - exp(-x)
    else {
      nn <- length(x)
      xx <- 1 + shape * x
      xxpos <- xx[xx > 0 | is.na(xx)]
      scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
      d <- numeric(nn)
      d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
        (1/shape + 1) * log(xxpos)
      d[xx <= 0 & !is.na(xx)] <- -(1e+6)
    }
    return(d)
  }
  
  start <- list()
  start$scale <- sqrt(6 * var(x))/pi
  start$loc <- mean(x) - 0.58 * start$scale
  p = ncol(z)

  zz = cbind(1,z)
  tvec = rep(0,p+3)
  tvec[1] = start$loc
  tvec[p+2] = start$scale
  return( optim(tvec, l2gev, lambda = lambda  , method = "BFGS")$par ) 
}
  
### simulation 

n = 100
true_theta=c(100,30,0.1)
p= 100
true_beta = rep(0,p)
# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3


set.seed(1)
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
z = matrix(rnorm(n*p), n, p)
x = z%*%true_beta + eps 
# ridge
gevreg(x,z, lambda = 10)
# function in gev package
fgev(x = x, nsloc = as.data.frame(z))$estimate

