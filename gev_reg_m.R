rm(list = ls())
library(evd)
# multiple GEV optim for the global trend
  # multiple sites
  # tvec( stationary: site1_para, site2_para, ..., 
#        nonstationary: common reg model)
  
  n = 100
  true_theta=c(100,30,0.1)
  p= 10
  true_beta = rep(0,p)
  # optim control
  optim_controlList = list()
  optim_controlList$maxit = 1e+3
  
  
  set.seed(1)
  xlist = zlist = list()
  ns = 10
  
  for (i in 1:ns)
  {
    z = matrix(rnorm(n*p), n, p)
    eps = rgev(n,loc=true_theta[1], 
               scale=true_theta[2], 
               shape=true_theta[3])
    xlist[[i]] = z%*%true_beta + eps 
    zlist[[i]] = z
  }
  

  gevreg_m = function(xlist, zlist, lambda = 0)
  {
    p = ncol(zlist[[1]])
    tvec = rep(0, ns*3 + p)
    ns = length(xlist)
    
    
    l2gev_m = function (tvec, lambda, xlist, zlist)
    {
      ns = length(xlist)
      v1 = 0
      for ( i in 1: length(xlist))
      {
        x = xlist[[i]]
        z = zlist[[i]]
        loc.vec.reg = drop(z%*%tail(tvec, p))
        loc.vec = tvec[3*(i-1)+1] + loc.vec.reg
        sc.vec = tvec[3*(i-1)+2]
        sh.vec = tvec[3*(i-1)+3]
        v1 = v1 - sum(lgev(x, loc = loc.vec, 
                           scale = sc.vec, shape = sh.vec))  
      }
      # loglikelihood
      
      # regularization
      v2 = lambda*sum(tail(tvec,p)^2)
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
    
    
    for ( i in 1:ns)
    {
      x = xlist[[i]]
      start <- list()
      start$scale <- sqrt(6 * var(x))/pi
      start$loc <- mean(x) - 0.58 * start$scale
      tvec[3*(i-1)+1] = start$loc
      tvec[3*(i-1)+2] = start$scale
    }
    
    
    return( optim(tvec, l2gev_m, lambda = lambda, 
                  method = "BFGS", 
                  control = optim_controlList,
                  xlist = xlist,
                  zlist = zlist)$par) 
  }
  
 fitvec = gevreg_m(xlist, zlist, lambda = 1)
 matrix(head(fitvec, 3*ns), ns, 3, byrow = T)
 tail(fitvec,p)
