# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.stats import genextreme

def gevreg_m(xlist, zlist, lambda =0 ):
    p = zlist[1].shape[1]
    ns = len(xlist)
    tvec = np.repeat(0, ns*3+p)
    
    def gev(x, loc = 0, scale = 1, shape = 0):
        if (scale <= 0) :
            return ( -1e+6)
        x = (x - loc)/scale
        if (shape = 0)
    

n = 100
true_beta = np.array([100,30,0.1], dtype = float)
p = 10
true_beta = np.repeat(0,p)
ns = 10
xlist = zlist = []



for i in range(ns):
    z = np.random.normal(size = n*p)
    z = z.reshape(n,p)
    x = genextreme.rvs(loc = true_beta[0],
               scale = true_beta[1],
               c = true_beta[2],
               size=n)
    xlist.append(x)
    zlist.append(z)

z.shape[1]
