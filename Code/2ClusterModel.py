#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 12:11:46 2017

@author: Jonathan Peters
"""

import numpy as np
import pyphi as phi
import matplotlib.pyplot as plt
import itertools
from numba import jit

states = np.array(list(itertools.product([0,1], repeat = 5)))


@jit
def transform_tpm(tpm):
    #this function changes a transition probability matrix from the 2^n x n format back to the 2^n x 2^n format.
    
    nrows = np.size(tpm,0);
    ncols = np.size(tpm,1);
    
    if nrows==ncols:
        print('do something');
    
    elif np.log2(nrows) == ncols:
        
        options = np.array(list(itertools.product([0,1], repeat = 5)))

        new = np.ones([nrows, nrows]);
        
        for i in range(nrows):
            for j in range(nrows):
                for k in range(ncols):
                    adj = (tpm[i,k] if options[j,k] else 1-tpm[i,k])
                    new[i,j] *= adj              
        
        return new
    
@jit
def stat_dist(tpm):
    eps = 1e-4;
    eig, vec = np.linalg.eig(tpm.T);
    
    eig = np.sort(eig);
    
    if eig[0]>-1+eps and eig[-2] < 1-eps:
        ans =vec.T[0];
        #if ans[0] <0:
        #    ans = np.multiply(-1, ans)
        #print(abs(ans))
        scale = np.sum(ans)
        return np.divide(ans,scale)
    else:
        print('A unique stationary distribution does not exist')
        print(eig)
        print(tpm)


#%%

def buildTPM(cm, probs, k):
    nnodes = 5;
    
    
    #strongly connected
    sc= np.array([[cm[i,j]*cm[j,i] for i in range(nnodes)] for j in range(nnodes)])
    
    CM = np.zeros([8,8], dtype = int)
    CM[:5,:5]=cm
    CM[5:, :5] = np.ones([3,5])
    
    states = np.array(list(itertools.product([0,1], repeat = 5)))
    projects = np.array(list(itertools.product([0,1], repeat = 3)))
    TPM = np.empty([2**8, 8], dtype = float);
    TPM[:,5] = np.multiply(probs[0],np.ones([2**8]))
    TPM[:,6] = np.multiply(probs[1],np.ones([2**8]))
    TPM[:,7] = np.multiply(probs[2],np.ones([2**8]))
    
    
    
    i = 0;
    
    
    
    for project_index in range(8):
        project = projects[project_index, ::-1]
        for row in range(32):
            
            s = states[row,::-1];
            ns = 1-s
            
            s_off = cm*ns;
            sc_off = sc*ns;
            sc_on = sc*s;

            #first determine how a project gets allocated
            if (project == [0,0,0]).all():
                A,B,C,D,E = 0,0,0,0,0
    
            elif (project == [1,0,0]).all():
                A = ns[0]
                B = 0;
                C = ~(s_off[2,0]|s_off[2,3])&ns[2]
                D = ~(s_off[3,0])&ns[3]
                E = 0;
                
            elif (project == [0,1,0]).all():
                A = 0;
                B = ns[1]
                C = ~(s_off[2,1]&s_off[2,4])&ns[2]
                D = 0;
                E = ~(s_off[4,1])&ns[4]
        
            elif (project == [1,1,0]).all():
                A = ~s_off[0,2]&ns[0]&(sc_off[0,1]|sc_off[0,4])
                B = ~s_off[1,2]&ns[1]&(sc_off[1,0] | ( ~(s_off[1,0]&s_off[1,4])&sc_off[1,3]) )
                C = ns[2];
                D = ~(s_off[3,2] | (s_off[3,0]&s_off[3,1]) | (s_off[3,4]&s_off[3,0]))&ns[3]&(sc_off[3,1] | sc_off[3,4])  
                E = ~(s_off[4,2] | (s_off[4,0]&s_off[4,1]))&ns[4]&( sc_off[4,0] | (~(s_off[4,1]&s_off[4,3])&sc_off[4,3]) )
                
                
            elif (project == [0,0,1]).all():
                A = 0;
                B = 0;
                C = 0;
                D = ~(s_off[3,4])&ns[3];
                E = ns[4];
                
            elif (project == [1,0,1]).all():
                A = ~(s_off[0,3])&sc_off[0,4]*ns[0];
                B = 0;
                C = ~(s_off[2,3] | (s_off[2,0]&s_off[2,4]) )&sc_off[2,4]*ns[2];
                D = ns[3];
                E = ~(s_off[4,3])&(s_off[4,0]|s_off[4,2])&ns[4];
                
            elif (project == [0,1,1]).all():
                A = 0;
                B = ~(s_off[1,4])&sc_off[1,3]&ns[1];
                C = ~( s_off[2,4] | (s_off[2,1]&s_off[2,3]))&sc_off[2,3]&ns[2];
                D = ~(s_off[3,4])&(s_off[3,1]|s_off[3,2])&ns[3];
                E = ns[4];
                
            elif (project == [1,1,1]).all():
                
                A = sc_off[0,4]&ns[0];
                B = ~(s_off[1,0]&s_off[1,4])&sc_off[1,3]&ns[1];
                C = ~( (s_off[2,0]&s_off[2,4]) | (s_off[2,1]&s_off[2,3]) | s_off[2,3]&s_off[2,4])&ns[2]&(sc_off[2,3]|sc_off[2,4])
                D = ~(s_off[3,0]&s_off[3,4])&ns[3]&(( sc_off[3,1] | sc_off[3,4] ) | ~(s_off[3,2]&s_off[3,4])&sc_off[3,2])
                E = (sc_off[4,0]&ns[4]) | ~(s_off[4,1]&s_off[4,3]) & ( sc_off[4,3] | sc_off[4,2])&ns[4]; 
            
            #next determine how projects are completed.
            
            #determine the chance of A completing its task
            if s[0]:
                A = 1-k[0];
                if sc_on[0,2]:
                    A-=k[2];
                if sc_on[0,3]:
                    A-=k[3];
            
            #determine the chance of B completing its task
            if s[1]:
                B = 1-k[1];
                if sc_on[1,2]:
                    B-=k[2];
                if sc_on[1,4]:
                    B-=k[4];
            
            #determine the chance of C completing its task
            if s[2]:
                C = 1-k[2];
                if s[0]&s[1]:
                    C-=(k[0]+k[1])
            
            #determine the chance of D completing its task
            if s[3]:
                D = 1-k[3];
                
            #determine the chance of E completing its task
            if s[4]:
                E = 1-k[4];
            
            
            
            TPM[row+32*project_index, 0:5] = A,B,C,D,E
            print(i, row+32*project_index,s, project, A,B,C,D,E, np.round(TPM[row+32*project_index, :5],3))
            i+=1
            
    
    
    print(TPM[100])
    return TPM,CM
    
                                                 
cm = np.zeros([5,5],dtype = int)
cm[0,1], cm[0,2], cm[1,0], cm[1,2], cm[2,0], cm[2,1] = 1,1,1,1,1,1
cm[3,4] = 1
cm[4,3]=1
k = [.2, .2, .1, .05, .05]


#%%
states = np.array(list(itertools.product([0,1], repeat = 8)))


TPM,CM = buildTPM(cm, (.5,.5,.5),k )




#%%
net = phi.Network(TPM, CM);
sub = phi.Subsystem(net, [0,1,1,0,1,1,0,0], range(8))
cmplxs = phi.compute.condensed(net, [0,1,1,0,1,1,0,0])
#%%

dat1 = np.zeros([2**8])
dat2 = np.zeros([2**8])




net = phi.Network(TPM, CM)
for state in states:
    state = state[::-1]
    dat = phi.compute.main_complex(net, state).phi
    tmp = phi.compute.condensed(net, state)
    
    print(state, dat, len(tmp))






