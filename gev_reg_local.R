rm(list = ls())
library(evd)
# multiple GEV optim for the local trend 
# this code is used for 2d spline of the GEV regression 
# multiple sites
  # tvec( stationary: site1_para (loc, loc_reg, scale, shape), 
  #                   site2_para (loc, loc_reg, scale, shape),
  # , ..., 
  #        nonstationary: common reg model)
  
  n = 100
  true_theta=c(100,30,0.1)
  p= 5
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

  gevreg_l = function(xlist, zlist, lambda = 0)
  {
    p = ncol(zlist[[1]])
    # para: (loc, loc_reg, scale, shape)
    tvec = rep(0, ns*(3 + p))
    ns = length(xlist)
    idx.vec = rep(F, (3+p)*ns)
    for ( i in 1: length(xlist))
    {
      idx = ((3+p)*(i-1)+2):((3+p)*(i-1)+1+p)
      idx.vec[idx] = T
    }
    #
    
    
    l2gev_m = function (tvec, lambda, xlist, zlist, idx.vec)
    {
      ns = length(xlist)
      v1 = 0
      for ( i in 1: length(xlist))
      {
        x = xlist[[i]]
        z = zlist[[i]]
        idx = ((3+p)*(i-1)+2):((3+p)*(i-1)+1+p)
        loc.vec.reg = drop(z%*%tvec[idx])
        loc.vec = tvec[(3+p)*(i-1)+1] + loc.vec.reg
        sc.vec = tvec[(3+p)*(i-1)+p+2]
        sh.vec = tvec[(3+p)*(i-1)+p+3]
        v1 = v1 - sum(lgev(x, loc = loc.vec, 
                           scale = sc.vec, shape = sh.vec))  
      }
      # loglikelihood
      (3+p)*(0:(ns-1))
      # regularization
      v2 = lambda*sum(tvec[idx.vec]^2)
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
      tvec[(3+p)*(i-1)+1] = start$loc
      tvec[(3+p)*(i-1)+1+p+1] = start$scale
    }
    
    
    return( optim(tvec, l2gev_m, lambda = lambda  , 
                  method = "BFGS", xlist = xlist,
                  zlist = zlist, idx.vec = idx.vec)$par) 
  }
  
  
  
 a = round(gevreg_l(xlist, zlist, lambda = 0.1),3)
 matrix(a,ncol = p+3, byrow = T)
  
