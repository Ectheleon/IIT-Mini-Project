#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 09:18:12 2017

@author: Jonathan Peters
"""

import numpy as np
import pyphi
import scipy.optimize as so
from numba import jit


@jit
def func(x,state =(0,0,0), n = 3):
    
    
    if len(x)==n*2**n:
        tpm = np.reshape(x, (2**n,n));
        net = pyphi.Network(tpm);
        subsystem = pyphi.Subsystem(net, state,range(net.size))
        phi = pyphi.compute.big_phi(subsystem);
        #print(how_close_to_int(x), phi)
        
        return phi
    else:
        print("Invalide input, not of the correct size to change into a TPM");
        
def how_close_to_int(x):
    return np.sum(np.array([xx*(1-xx) for xx in x]));

#%%

func([0,1,.5,.5,0.1,0.9,1,0],(0,0), n=2)
#%%
n = 3
x0= np.array([0.5 for i in range(n*2**n)])

bounds = np.array([[0,1] for i in range(n*2**n)])

res = so.differential_evolution(func, bounds, maxiter = 10, tol = 0.1, disp = True)

#%%
def dist_from_bounds(x):
    return np.array([np.minimum(0.999-xx, xx-0.001) for xx in x])

n = 3;
state = (0,0,0)

ntests = 500
per_size = 0.2
start = np.array([0.5 for i in range(n*2**n)])
per_size = dist_from_bounds(start);


phis = []
vec = []
vec.append(start)
dat = []

x = start;
fold = func(x, state);
phis.append(fold);
dat.append(how_close_to_int(x));

for i in range(ntests):
    if np.min(dist_from_bounds(x)- per_size)<0:
        per_size = dist_from_bounds(x);
    perturbation = np.random.uniform(-1, 1, n*2**n)
    can = x +per_size*perturbation/2 #candidate
    
    fcan = func(can, state);
    
    if fold ==0:
        aprob = 1;
    else:
        if fcan<fold:
            aprob = fcan/(4*fold)
        else:
            aprob = 1;
    
    u = np.random.uniform(0,1)
    if u < aprob:
        print("{0:0.4}".format(how_close_to_int(can)), "{0:0.4}".format(fold))
        x = can
        fold = fcan;
        phis.append(fold);
        vec.append(x)
        dat.append(how_close_to_int(can))

for i in np.reshape(vec[-1], (2**n, n)):
    print (i)







