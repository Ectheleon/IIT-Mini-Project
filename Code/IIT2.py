#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:23:45 2017

@author: jonathan
"""
import numpy as np
#%%

def convert_to_subset( num , subset ,base=2):
    '''
    An integer describes the state of a number of nodes, the state of the ith node
    being defined as the ith digit of the integer when in base 'base'. 
    
    When talking about a subset of nodes, its convenient to condense this notation.
    Thus the ith digit now refers to the ith node in the subset. 
    
    This function converts a integer referring to the state of the system to the 
    reduced state which only refers to the subset.
    
    For example: 5 nodes, 0, 1, 2, 3, 4.
    
    The current state is 01101 which is described by the integer 13.
    
    For the subset 0, 2,3 the corresponding state is 111, and thus the integer 
    is 7. This function changes 13 to 7 given the subset [0, 2, 3]
    '''
    
    subset_nodes = np.array(sorted(subset), dtype = int);
    subset_size = np.size(subset_nodes);
    
    return np.sum(bitget(num,subset_nodes[j])*base**j for j in xrange(subset_size));
    

def convert_from_subset( num, subset, base=2):
    '''
    This function is the inverse of the previous one. Given an integer which
    describes the state of a subset, this will give an integer corresponding
    to an overall state, where all elements of the subset have their correct state.
    
    This answer is not unique. Given the subset [0,2,3] in [0,1,2,3,4], the binary
    string of 111 could equivelantly be translated to 11111 or 10110. 
    
    For convenience later, we select the integer where the state of all nodes not
    within the subset is 0.
    '''
    
    subset_nodes = np.array(sorted(subset), dtype = int);
    subset_size = np.size(subset_nodes);
    
    return np.sum(bitget(num, j)* base**subset_nodes[j] for j in xrange(subset_size));    
#%%
    
    

def bitget(integer, i, base=2):
    '''This function allows integers to be viewed as an array of base n numbers.
    Bitget will return the value of the i'th digit of integer when expressed in 
    base: 'base'. By default we assume to be working with binary numbers.
    '''
    if base==2:
        return (integer>>i)&1
    #The 'else' case we will neglect for now.
        


def perturb( subset , base = 2):
    '''subset is a set of integers, each corresponding to a node. 
        
    By default, it is assumed that each node has two states, hence 'base' = 2.
    
    With n nodes, each haveing 'base' states, it follows that there are base^n
    total states.
    
    perturb returns a vector of length base^n, and an integer. The array contains
    the locations the base^n states which can all occur with equal liklihood. The
    integer contains the probability of them happening.
    '''
    
    nodes = list(sorted(subset));
    nnodes = len(subset);
        
    locations = np.array([np.sum([bitget(i,j)*base**nodes[j] for j in xrange(nnodes)]) for i in xrange(base**nnodes)],dtype = int)
    
    return locations, 1.0/(base**nnodes);
        
def uncon_cause_repertoire( n, base = 2 ):
    '''
    This function returns the unconstrained cause repertoire for a system of n nodes;   
    This will be the uniform distribution across all possible states.
    '''
    return np.repeat(1.0/n, n);

def cause_repertoire( current, past, state, tpm, base=2):
    '''
    This function caluclates the distribution of possible past states given 
    particular information of the present state. 
    
    current: is a set of integers referring to the nodes we know something about.
    past: a list of integers referring to the nodes we want to find out about.
    state: is the current state which the nodes defined by 'current' are in.
    tpm: is the 2^n x n transition probability matrix.
    '''
    
    #With respect to the tpm, the only columns we are interested in are those 
    #which are indexed by the label of the set of current nodes.
    cols = np.array(sorted(current), dtype = int);
    
    #With respect to the tpm, the only rows we are interested in are those which
    #correspond to possible 
    
def condense_tpm( denominator , tpm , base = 2 ):
    nodes = np.array(sorted(denominator), dtype = int);
    nnodes = len(nodes);
    
    condensed_tpm = np.zeros(base**nnodes, np.size(tpm, 1))
    
    for i in xrange(base**nnodes):
        
        
    
    
    
    

tmp = set([2,1]);
perturb(tmp)