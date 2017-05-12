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

target_partition = set([0,1]);
known_partition = set([1,2,3]);

known_state = np.array([1,0,0]);


#%% Trial the function

partial_probability(target_partition, known_partition, known_state, tpm);


#%%

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
    
    whole_set = set(range(nnodes));
    
    target_size = len(target_partition);
    target_compliment = whole_set.difference(target_partition);
    target_compl_nnodes = len(target_compliment);
    target_compl_nstates = target_compl_nnodes;
    target_list = list(sorted(target_partition));
    
    known_size = len(known_partition);
    known_compliment = whole_set.difference(known_parition);
    known_compl_nnodes = len(known_compliment);
    known_compl_nstates = 2**known_compl_nnodes;
    known_list = list(sorted(known_partition));
    
    
    if target_compl_nnodes ==0:
        perturbation = 0;
    else:
        perturbation = np.zeros(known_compl_nstates)
        
        for j in xrange(known_comp_nstates):
            perturbation[i] = 
        
    
    
    
    