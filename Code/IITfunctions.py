#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 15:46:28 2017

@author: Jonathan Peters
"""

import IIT2 as iit
import numpy as np
from emd import emd

#%%
reload(iit)

def cei( subset , state , tpm , base = 2 ):
    '''
    This function calculates the cause-effect imformation resulting
    from the state of a particular subset being known. 
    
    subset is the set of nodes for which the state is known.
    state is an integer which describes the state of these nodes.
    tpm is the 2^n x n transition probability matrix.
    '''
    
    nnodes = np.size(tpm, 1);
    full_set = set(range(nnodes));
    
    #print(nnodes);
    
    f_uncon = iit.uncon_effect_repertoire(tpm, base);
    p_uncon = iit.uncon_cause_repertoire(nnodes, base);
    
    
    f = iit.effect_repertoire(subset, full_set, state, tpm, base);
    p = iit.cause_repertoire(subset, full_set, state, tpm, base);
    
    '''
    cause_information = iit.EMD1(p_uncon, p);
    effect_information = iit.EMD1(f_uncon, f);
    
    print('ci',cause_information)
    print('ei',effect_information)
    
    cause_information = iit.EMD2(p_uncon, p);
    effect_information = iit.EMD2(f_uncon, f);

    print('ci',cause_information)
    print('ei',effect_information)
    
    Dist = np.array([[iit.hamming(i,j) for i in xrange(2**nnodes)] for j in xrange(2**nnodes)], dtype=int)  
    print (Dist)
    '''
    
    d = np.array(range(2**nnodes))
    locs = ((d[:,None] & (1 << np.arange(nnodes-1, -1, -1))) > 0).astype(int)

    cause_information = emd(locs, locs, p_uncon, p, distance='cityblock');
    effect_information = emd(locs, locs,f_uncon, f, distance = 'cityblock');

    
    return np.minimum(cause_information, effect_information);

def pre_phi_past(current, past, current_part, past_part, state, tpm, base=2):
    '''
    This is a preliminary function which makes the calculation of small_phi 
    cleaner, although it itself has no formal name.
    
    Given a bipartition of the current set and the cause set, this function
    calculates the distance between the repertoire specified by the whole 
    mechanism, and the repertoire specified by the partitioned mechanism
    
    current is the set of nodes which we know something about. 
    state is what we know about those nodes
    current_part is a subset of current. The other subset is its compliment
    
    past is the set of nodes from the past we are interested in.
    past_part is a subset of past, the other subset is its compliment.
    '''
    
    whole_rep = iit.cause_repertoire(current, past, state, tpm, base);
    
    part1_state = iit.convert_to_subset(state, current_part, base);
    part1_rep = iit.cause_repertoire(current_part, past_part, part1_state, tpm,base);
    
    part2_state =iit.convert_to_subset(state, current - current_part, base);
    part2_rep = iit.cause_repertoire(current - current_part, past - past_part, part2_state, tpm,base);
    
   
    partitioned_rep = iit.multiply_repertoires(past_part, past-past_part, part1_rep, part2_rep, base);
        
    past_nnodes = len(past);
    d = np.array(range(2**past_nnodes))
    locs = ((d[:,None] & (1 << np.arange(past_nnodes-1, -1, -1))) > 0).astype(int)
    
    return emd(locs, locs, whole_rep, partitioned_rep, distance='cityblock');

def pre_phi_future(current, future, current_part, future_part, state, tpm, base=2):
    whole_rep = iit.effect_repertoire(current, future, state, tpm, base);
    
    part1_state = iit.convert_to_subset(state, current_part, base);
    part1_rep = iit.effect_repertoire(current_part, future_part, part1_state, tpm,base);
    
    part2_state =iit.convert_to_subset(state, current - current_part, base);
    part2_rep = iit.effect_repertoire(current - current_part, future - future_part, part2_state, tpm,base);
    
   
    partitioned_rep = iit.multiply_repertoires(future_part, future-future_part, part1_rep, part2_rep, base);
        
    future_nnodes = len(future);
    d = np.array(range(2**future_nnodes))
    locs = ((d[:,None] & (1 << np.arange(future_nnodes-1, -1, -1))) > 0).astype(int)
    
    return emd(locs, locs, whole_rep, partitioned_rep, distance='cityblock');

def phi_MIP_cause(current, past, state, tpm, base = 2):
    
    min_phi = np.inf;
    MIP = set()
    size1 = len(current)
    size2 = len(past);
    count = 0;
    for i in iit.powerset(past):
        for j in iit.powerset(current):
            if count >0 and count < base**(size1+size2)-1:
                temp = pre_phi_past(current, past, j, i, state, tpm, base);
                
                if temp < min_phi:
                    min_phi = temp;
                    MIP = [j,i];
            count+=1;
                
    if min_phi<np.inf:
        return min_phi, MIP        
    else:
        return 0, MIP

def phi_MIP_effect(current,future, state, tpm, base = 2):
    min_phi = np.inf;
    MIP = set()
    
    size1 = len(current)
    size2 = len(future);    
    count = 0;
    for i in iit.powerset(future):
        for j in iit.powerset(current):
            if count >0 and count < base**(size1+size2)-1:
                temp = pre_phi_future(current, future, j, i, state, tpm, base);
                
                if temp < min_phi:
                    min_phi = temp;
                    MIP = [j,i];
            count+=1;
                
    if min_phi<np.inf:
        return min_phi, MIP      
    else: 
        return 0, MIP

def phi_MIP(current, past, future, state, tpm, base=2):
    '''
    This function determines the Minimum Information Partition given the 
    state: 'state' on the set of nodes 'current' projected into both the past
    and the future. 
    
    The value of phi^MIP is returned as well as the MIP, in fractional form.
    '''
    
    a,b = phi_MIP_cause(current, past, state, tpm, base)
    c,d = phi_MIP_effect(current, future, state, tpm, base)
    
    return np.minimum(a,c), [b, d]

def phi_MAX_cause(current, full_set, state, tpm, base = 2):
    '''
    calculates the core cause set of the set 'current' being is state 'state'.
    Returns the value of phi^max_cause as well as the core cause.
    '''
    phi_max = 0;
    Core_Cause = set()
    
    for i in iit.powerset(full_set):
        temp, candidate = phi_MIP_cause(current, i, state, tpm, base)
        if temp > phi_max:
            phi_max = temp;
            Core_Cause = i;

    return phi_max, Core_Cause;

def phi_MAX_effect(current, full_set, state, tpm, base = 2):
    '''
    calculates the core effec set of the set 'current' being is state 'state'.
    Returns the value of phi^effect as well as the core effect.
    '''
    phi_max = 0;
    Core_Effect = set()
    
    for i in iit.powerset(full_set):
        temp, candidate = phi_MIP_effect(current, i, state, tpm, base)
        if temp > phi_max:
            phi_max = temp;
            Core_Effect = i;

    return phi_max, Core_Effect;
    
def phi_MAX(current, full_set, state, tpm, base=2):
    '''
    determines the Maximally Irreducible Cause-Effect Repertoire, and the 
    corresponding value of phi_MAX.
    '''
    
    a,b = phi_MAX_cause(current, full_set, state, tpm, base=2);
    c,d = phi_MAX_effect(current, full_set, state, tpm, base=2);
    
    return np.minimum(a,c), [b,d]
    
    
    
#%%
tpm = np.array([[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0], [1,1,1], [1,0,1],[1,1,0]]);


past = set([0,1,2])
current = set([0,1,2])
future = set([0,1,2])

state = 1


phi_MAX(set([0]), current, 1,tpm)
#phi_MIP( current,past, future, state ,tpm)


#%%
current = set([1,2])
past = set([0,2])

current_part = set([2])
past_part = set([])

state = 0

pre_phi(current, past, current_part, past_part, state, tpm)



#%%
subset = set([0]);
state = 1;
tpm = np.array([[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0], [1,1,1], [1,0,1],[1,1,0]]);


iit.uncon_effect_repertoire(tpm)

ans = cei(subset, state, tpm);

print(ans);

#%%

