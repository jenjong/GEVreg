# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import genextreme

n = 3
ns = 2
p = 2
true_theta = np.array([100,30,0.1], dtype = float)
true_beta = np.repeat(0.0,p)
xlist = []
zlist = []
for i in range(ns):
    z = np.random.normal(size = n*p)
    z = z.reshape(n,p)
    x = genextreme.rvs(loc = true_theta[0],
               scale = true_theta[1],
               c = true_theta[2],
               size=n)
    xlist.append(x)
    zlist.append(z)

def gevreg_m(xlist, zlist, lambda =0 ):
    p = zlist[1].shape[1]
    ns = len(xlist)
    tvec = np.repeat(0.0, ns*3+p)
    

    def lgev(x, loc = 0, scale = 1, shape = 0):
        if (scale <= 0) :
            return ( -1e+6)
        x = (x - loc)/scale
        if (shape == 0):
            d = np.log(1/scale) -x -np.exp(-x)
        else:
            nn = len(x)
            xx = 1 + shape*x
            xxpos = xx[xx >0 | np.isnan(x)]
            scale = np.repeat(scale, nn)[ (xx>0) |(np.isnan(x))]
            d =  np.repeat(0.0, nn)
            d[xx > 0 | np.isnan(xx)] =  np.log(1/scale)  \
            - np.power(xxpos,(-1/shape)) \
            -(1/shape + 1) * np.log(xxpos)
            d[(xx <= 0) & (~np.isnan(xx))] = -(1e+6)
        return(d)
    

    x = xlist[1]        
    loc = 100.0
    scale = 20.0
    shape = 0.1    
    lgev(xlist[1], loc = loc, scale = scale, shape = shape)
    lam = 1.0
    # 



    for i in range(ns):
        x = xlist[i]
        start_scale = np.sqrt(6.0 * np.var(x))/np.pi
        start_loc = np.mean(x) - 0.58 * start_scale
        tvec[3*i] = start_loc
        tvec[3*i+1] = start_scale
    
    
    args = {'lam': lam, 'xlist': xlist, 'zlist':zlist}
    
    
    def l2gev_m (tvec, args):
        lam, xlist, zlist = args['lam'], args['xlist'], args['zlist']
        ns = len(xlist)
        v1 = 0
        for i in range(ns):
            x = xlist[i]
            z = zlist[i]
            loc_vec_reg = np.matmul(z,tvec[-p:(ns*3+p)])
            loc_vec = tvec[3*i] + loc_vec_reg
            sc_vec = tvec[3*i+1]
            sh_vec = tvec[3*i+2]
            v1 = v1 - np.sum(lgev(x, loc = loc_vec, scale = sc_vec,
                           shape = sh_vec))
        v2 =  lam*np.sum(np.power(tvec[-p:],2))
        v = v1 + v2
        return (v)
    
    l2gev_m (tvec, args)
  
    fit = minimize(l2gev_m, x0 = tvec, method='BFGS',args = args )

