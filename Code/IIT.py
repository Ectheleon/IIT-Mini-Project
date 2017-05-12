#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:08:02 2017

@author: Jonathan Peters
"""

import numpy as np
import scipy as sp

#%% Preliminaries

filepath = 'External Software/iit/'

temp = sp.io.loadmat(filepath+'tpm_matrix_M1_K0.10_0.10_0.10_P0.50_0.02.mat');
TPM = temp['tpm'];
tpm = convert_tpm(TPM);

target_partition = set([1]);
known_partition = set([1]);

known_state = np.array([1,0,0]);




#%%

def bitget(integer, bit):
    return (integer>>bit)&1;


def convert_tpm(tpm):
    nstates = np.size(tpm,0)
    nnodes = int(np.log2(nstates));
    
    new_tpm = np.zeros([nstates, nnodes]);
    
    for i in xrange(nstates):
        for j in xrange(nnodes):
            for k in xrange(nstates):
                if k& (1<<j):
                    new_tpm[i,j] += tpm[i,k];
    return new_tpm;


def partial_probability(target_partition, known_partition, known_state, TPM):
    '''
    This function caluclates the probability distribution of states which the
    a partition of nodes (of_partition) might take at a timestep, given that 
    at the previous time step we know that a distinct partition (given_partition)
    has the state 'known_state'.
    
    The size of the system should be contained within the matrix: TPM which 
    refers to the transition probability matrix. This matrix should be of 
    dimensions: 2^n x n, where n is the number of nodes, and hence 2^n is the
    number of distinct states that the network my be in. (Each node has 2 states)
    
    target_partiton: should be a python set object, which should contain the indices 
    of the nodes which are in the partition of question.
    
    known_partition: should similarly be a python set object as before.
    
    known_state should be a vector of length n. This vector should contain only
    1's or 0's and indicates whether the corresponding node is turned on or off.
    '''
    
    nnodes = np.size(TPM, 1);
    nstates = np.size(TPM, 0);
    
    pows_2 = np.array([2**i for i in range(nnodes)]);
    
    whole_set = set(range(nnodes));
    
    target_size = len(target_partition);
    target_list = list(sorted(target_partition));
    target_compliment = whole_set.difference(target_partition);
    target_compl_nnodes = len(target_compliment);
    target_compl_nstates = 2**target_compl_nnodes;
    target_compl_list = list(sorted(target_compliment));
    
    known_size = len(known_partition);
    known_list = list(sorted(known_partition));
    known_compliment = whole_set.difference(known_partition);
    known_compl_nnodes = len(known_compliment);
    known_compl_nstates = 2**known_compl_nnodes;
    known_compl_list = list(sorted(known_compliment));
    

    if known_compl_nnodes ==0:
        known_perturbations = np.zeros(0);
    else:
        #temp = np.array([ 2**i for i in known_compl_list], dtype = int);
        known_perturbations = np.array([np.sum([pows_2[known_compl_list[i]] * bitget(j,i) for i in range(known_compl_nnodes)]) for j in range(known_compl_nstates)])
    
    if target_compl_nnodes ==0:
        target_perturbations = np.zeros(0);
    else:
        #temp = np.array([ 2**i for i in known_compl_list], dtype = int);
        target_perturbations = np.array([np.sum([pows_2[target_compl_list[i]] * bitget(j,i) for i in range(target_compl_nnodes)]) for j in range(target_compl_nstates)])
        
    
    #calculate the 'prob' vector.
    prob = np.zeros([2**target_size]);
    for i in xrange(2**target_size):
        prob[i] = 1;
        var1 = np.sum([pows_2[target_compl_list[j]] * bitget(i,j) for j in range(target_compl_nnodes)]) 
        
        for k in xrange(known_size):
            var2 = known_state[known_list[k]];
            p_k = 0;
            
            for j in xrange(target_compl_nstates):
                var3 = target_perturbations[j];
                var4 = 1 + var1 + var3;
                
                if var2 ==1:
                    p_k+= TPM[var4, known_partition[k]];
                else:
                    p_k += (1-TPM[var4, known_partition[k]]);
            prob[i]+=p_k;
    
        
    
        #temp2 = np.array([bitget()])
        
        #temp2 = np.array([ ((2**i if bitget(i,j) else 1) for i in known_list) for j in range(known_size)], dtype = int)
        
    #return temp2
    '''  
    for i in xrange(known_comp_nstates):
        perturbations[i] = np.sum([temp[j] if bitget(j,i) for j in known_])
    
    
    
    perturbations = np.array([((2**i for i in known_list) if else 1) for j in known_compl_nstates], dtype = int)
    
    for j in xrange(known_comp_nstates):
        perturbations[i] = np.sum([2**i for i in known_list) ])'''
    
#%% Trial the function

ex = partial_probability(target_partition, known_partition, known_state, tpm);
    
    
    