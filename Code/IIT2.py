#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:23:45 2017

@author: jonathan
"""
import numpy as np
import scipy.optimize as so
import sys
import itertools 

def powerset(ofset):
    #returns the powerset of inputted set
    return [set(z) for z in itertools.chain.from_iterable(itertools.combinations(ofset, r) for r in range(len(ofset)+1))]

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
    
def hamming(int1, int2, base = 2):
    '''
    calulates the hamming distance between two integers when interpreting these
    integers as arrase of base 'base' numbers
    '''
    if base==2:
        diff =  (int1^int2);
        
        count = 0;
        while diff>0:
            if diff&1:
                count+=1;
            diff>>=1
    return count;
        
    #This code is only equipped for base 2 at the moment

def bitget(integer, i, base=2):
    '''This function allows integers to be viewed as an array of base n numbers.
    Bitget will return the value of the i'th digit of integer when expressed in 
    base: 'base'. By default we assume to be working with binary numbers.
    '''
    if base==2:
        return (integer>>i)&1
    #This code is only equipped for base 2 at the moment
        


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
    return np.repeat(1.0/(2**n), 2**n);

def uncon_effect_repertoire( tpm , base = 2 ):
    '''
    This function returns the unconstrained future repertoire for a system of 
    n nodes
    '''
    nnodes = np.size(tpm, 1);
    nstates = base**nnodes;
    
    #Store in the vector probs the unconstrained probability that the ith node
    #will be in state 1 after one step into the future.
    probs = np.array([np.divide(np.sum(tpm[:,i]), 1.0*nstates) for i in xrange(nnodes)]);
    
    #Use probs to calculate the unconstrained future repertoire.
    return np.array([np.prod([probs[i] if bitget(j,i) else 1-probs[i] for i in xrange(nnodes)]) for j in xrange(nstates)])



def cause_repertoire( current, past, state, tpm, base=2):
    '''
    This function caluclates the distribution of possible past states given 
    particular information of the present state. 
    
    current: is a set of integers referring to the nodes we know something about.
    past: a list of integers referring to the nodes we want to find out about.
    state: is the current state which the nodes defined by 'current' are in.
    This is the integer applying to the subset, not the full system.
    tpm: is the 2^n x n transition probability matrix.
    '''
    
    #With respect to the tpm, the only columns we are interested in are those 
    #which are indexed by the label of the set of current nodes.
    cols = np.array(sorted(current), dtype = int);
    
    #current_subset_nodes = np.array(sorted(current),dtype = int)

    full_set_size = np.size(tpm,1);

    con_tpm = condense_tpm(set(range(full_set_size)), past, tpm, base);

    vec = np.array([np.prod([i[cols[j]] if bitget(state, j) else 1-i[cols[j]] for j in xrange(len(cols))]) for i in con_tpm]);

    if sum(vec)!=0:
        return np.divide(vec, sum(vec));
    
        
def effect_repertoire( current , future , state , tpm , base = 2 ):
    '''
    This funciton calculates the future repertoire of a subset of nodes given
    a known state in a possibly different subset of nodes
    
    current: the set of nodes about which something is currently known
    future: the set of nodes for which we want to know the distribution of states
    for a future time step
    state: the integer corresponding to the state of the nodes in 'current'
    tpm: the 2^n x n transition probability matrix.
    '''
        
    future_nnodes = len(future);
    future_nstates = base**future_nnodes;
    future_nodes = np.array(sorted(future),dtype = int)
    
    total_nnodes = np.size(tpm,1);
    
    #Store in the vector probs the unconstrained probability that the ith node
    #will be in state 1 after one step into the future.
    con_tpm = condense_tpm(set(range(total_nnodes)), current, tpm, base);
    
    probs = con_tpm[state];
    
    return np.array([np.prod([probs[future_nodes[i]] if bitget(j,i) else 1-probs[future_nodes[i]] for i in xrange(future_nnodes)]) for j in xrange(future_nstates)])
    
    
    
def condense_tpm( full_set, denominator , tpm , base = 2 ):
    '''vec
    When looking at only subsets of the full system, we often do not need the entire
    transition probability matrix. This function condensed the tpm into the 
    appropriate form for a subset of nodes.
    
    '''
    np.set_printoptions(threshold=np.inf)

    nodes = np.array(sorted(denominator), dtype = int);
    nnodes = len(nodes);
    
    
    total_nnodes = len(full_set);
    compliment = full_set - denominator;
    
    condensed_tpm = np.zeros([base**nnodes, total_nnodes])
    
    for i in xrange(base**nnodes):
        temp = convert_from_subset(i,denominator);
        
        perturbations , scaling = perturb(compliment);
        
        temp2 = np.add(perturbations, temp);
        
        condensed_tpm[i] = np.multiply(np.sum((tpm[j]) for j in temp2), scaling);
        
    return condensed_tpm;

#%% This function is unused, as it seems to have an error in it. It does not 
# return values of the EMD consistent with results from Tononi's paper. Moreover
# an external implementation of the EMD does agree. Hence the external softare
# is used.    
     
def EMD1(distr1, distr2, base=2):
    '''
    This function calculates the earth mover distance between two distributions
    of states which correspond to the same set of nodes.
    
    Given n nodes, there should be base**n different states, and this should be 
    the length of both input vectors. 
    
    The distance between each state is taken to be the hamming distance.
    '''
    
    #this should always apply
    if len(distr1) == len(distr2):
        nstates = len(distr1);
        #nnodes = int(np.log2(nstates)/np.log2(base));
        
        #define the distance matrix
        D = np.array([hamming(i,j,base) for i in range(nstates) for j in xrange(nstates)]);
        D = D.flatten();
        
        constraint = np.zeros([2*nstates, nstates**2])
        for i in xrange(nstates): 
            constraint[i, nstates*i:nstates*(i+1)] = np.ones(nstates);
            constraint[nstates:2*nstates, nstates*i:nstates*(i+1)] = np.eye(nstates);
        
        eq =np.concatenate([distr1, distr2]);
        
        ans = so.linprog(D, A_eq = constraint, b_eq = eq);
        
        #This should always be the case
        if ans.success:
            return ans.fun
        else:
            sys.exit('The earth mover distance failed to compute. This shouldn\'t happen!')

def multiply_repertoires(set1, set2, distr1, distr2,base=2):
    '''
    Suppose distr1 (which appllies to set 1) is a distribution over the states
    of n nodes, and distr2 applies to m nodes where m+n = N = total number of 
    nodes in the system. This function multiplies these distributions together
    and returns a distribution for the entire system.
    '''    
    
    total = set1|set2
    
    set1_nodes = np.array(sorted(set1),dtype = int)
    set2_nodes = np.array(sorted(set2),dtype = int)
    
    set1_nnodes = np.size(set1_nodes);
    set2_nnodes = np.size(set2_nodes);
    
    total_nnodes = set1_nnodes + set2_nnodes;
    
    temp = np.array([[convert_to_subset(convert_from_subset(i, set1,base)+convert_from_subset(j, set2, base),total) for i in xrange(2**set1_nnodes)] for j in xrange(2**set2_nnodes)])

    new_distr = np.ones(2**total_nnodes);
    for i in xrange(2**set1_nnodes):
        for j in xrange(2**set2_nnodes):
            new_distr[temp[j,i]] = distr1[i]*distr2[j]
    
    return new_distr
    
    

    
#%%

'''
tpm = np.array([[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0], [1,1,1], [1,0,1],[1,1,0]]);

a = cause_repertoire(set([0,1]), set([2]), 1, tpm)
b = cause_repertoire(set([2]), set([0,1]), 0, tpm)

aset = set([2])
bset = set([0,1])

multiply_repertoires(aset,bset,a,b)
'''
     
#%%
'''
tpm = np.array([[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0], [1,1,1], [1,0,1],[1,1,0]]);

print cause_repertoire(set([0]), set([0,1,2]), 1, tpm)

print cause_repertoire(set([1]), set([0,1,2]), 1, tpm)

print cause_repertoire(set([2]), set([0,1,2]), 1, tpm)

'''

#%%

'''
tpm = np.array([[0,0,0],[1,0,0],[1,0,1],[1,0,1],[0,0,1], [1,1,1], [1,0,0],[1,1,0]]);

distr1 = uncon_cause_repertoire(3);
distr2 = uncon_effect_repertoire(tpm);
print(distr1);
print(distr2);

print(EMD1(distr1, distr2));
'''
#%%
'''
tpm2 = condense_tpm( set([0,1,2]), set([0,2]), tpm )       
print(tpm2)
    
    
    

tmp = set([2,1]);
perturb(tmp)
'''