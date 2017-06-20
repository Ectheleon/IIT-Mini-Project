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
import scipy as sp;

states = np.array(list(itertools.product([0,1], repeat = 8)))
costs = np.array([.45,.3, .25]);

probs = np.array([0.5,.5,.5])

team_skill_map = np.zeros([5, 3])
team_skill_map[0,0] = 1;
team_skill_map[1,1] = 1;
team_skill_map[2,0], team_skill_map[2,1] = 1,1;
team_skill_map[3,2], team_skill_map[3,0] = 1,1;
team_skill_map[4,2], team_skill_map[4,1] = 1,1;

                                         
cm = np.zeros([5,5],dtype = int)
#cm[0,1] = 1
cm[0,2], cm[1,0], cm[1,2], cm[2,0], cm[2,1] = 1,1,1,1,1
cm[3,4] = 1
cm[4,3]=1
k = [.2, .2, .1, .05, .05]
cm[0,0] , cm[1,1], cm[2,2], cm[3,3],cm[4,4] = 1,1,1,1,1
cm[3,0],cm[3,1], cm[4,0], cm[4,1]=1,1,1,1

k = np.divide(np.ones([8]),8)


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

def connection_options(nedges):
    cm = np.eye(5);
    
    if (nedges<5) | (nedges>25):
        print('error, the network must have between 5 and 25 edges inclusive')
    
    nopts = sp.misc.comb(20,nedges)
    opts = np.empty([int(nopts), 5,5], dtype = int)

    #create a length 20 array of 0's;
    arr = np.zeros([20])
    #loop through all places to place ones.
    
    k = 0
    for i in itertools.combinations(range(20),nedges):
        arr.fill(0)
        #fill in the 1's
        for j in i:
            arr[j] = 1;
        #reshape
        temp = np.reshape(arr, [5,4]);
        
        cm[0,1:] = temp[0];
        cm[1,0] = temp[1,0];
        cm[1,2:] = temp[1,1:]
        cm[2,0:2] = temp[2,0:2]
        cm[2,3:] = temp[2,2:]
        cm[3,:3] = temp[3, :3]
        cm[3,4] = temp[3,3]
        cm[4,:4] =temp[4,:4]


        opts[k]= cm;
        k+=1
    return opts



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
                C = ~(s_off[2,1]|s_off[2,4])&ns[2]
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
            TPM = np.round(TPM, 3)
            #print(i,s, project, np.round(TPM[row+32*project_index, :5],3))
            i+=1
            
    
    
    return TPM,CM
    


def classical(TPM,skill_map, state, costs):
    #computes the cost of waste when the system transitions out of 'state'.
    waste = 0;
    
    old = state;
    ind = np.sum([2**i * old[i] for i in range(8)]);

    new = np.ceil(TPM[ind]);
    new_workers = new[:5]-old[:5];
    
    used_skills = new_workers.dot(skill_map);
    needed_skills = old[5:]
    
    
    #How many times are we doing this project?
    
    if np.sum(needed_skills)>0:
        if(np.sum(used_skills)>0):
            n = 0;
            
            while (used_skills-n*needed_skills>=0).all():
                n+=1;
            #now we know the project is being done n times.
            
            ok_excess = used_skills-(n-1)*needed_skills;
            if n>1:
                bad_excess = (n-2)*needed_skills;
            else:
                bad_excess = 0*needed_skills
            
            #we've allocated the project well, now the only waste is excess skills.
            for i in range(3):
                #we only count 10% of ok excess
                waste+= costs[i]*ok_excess[i]*0.1
                
                #we count the full amount of bad excess
                waste +- costs[i]*bad_excess[i]
                
        else:
            #now the question is, does there exist a team which should have 
            #taken up the new project?
            
            #do the required skills exist?
            if ((1-old[:5]).dot(skill_map)-needed_skills>=0).all():
                for i in range(3):
                    waste+=10*costs[i]*needed_skills[i];

    
    return waste;

def quantify_network(cm, probs, k):
    TPM, CM = buildTPM(cm,probs, k);
    net = phi.Network(TPM, CM)
    
    dat = np.zeros([2**8, 2])

    j= 0;
    for state in states:
        state = state[::-1]
        waste = (classical(TPM, team_skill_map, state,costs))
        cmplxs = phi.compute.condensed(net,state);
        
        phis = []
        for i in range(len(cmplxs)):
            phis.append(cmplxs[i].phi)
        dat[j,:] = np.array([waste, np.sum(phis)])
        #print(j,waste, phis)
        j+=1
        
    dist = np.real(stat_dist(transform_tpm(TPM)));
    
    return dist.dot(dat)
    
#%%
nedges = 10;
cms = connection_options(nedges);

dat = [];
for i in range(10):
    print('\n\nStarting calculation for new connectivity matrix')
    a,b = quantify_network(cms[i], probs, k)
    dat.append([a,b])
    print(cms[i], a, b)












#%%
TPM, CM = buildTPM(cm, probs, k)
net = phi.Network(TPM, CM);

dat = np.zeros([2**8, 2])

j= 0;
for state in states:
    state = state[::-1]
    waste = (classical(TPM, team_skill_map, state,costs))
    cmplxs = phi.compute.condensed(net,state);
    
    phis = []
    for i in range(len(cmplxs)):
        phis.append(cmplxs[i].phi)
    dat[j,:] = np.array([waste, np.sum(phis)])
    print(waste, phis)
    j+=1
    
#%%
plt.scatter(dat[:,0], dat[:,1])


#%%
full_TPM = transform_tpm(TPM)
avg_dist = np.real(stat_dist(full_TPM))