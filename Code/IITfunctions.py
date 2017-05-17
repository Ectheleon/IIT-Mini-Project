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



#%%

subset = set([0]);
state = 1;
tpm = np.array([[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0], [1,1,1], [1,0,1],[1,1,0]]);


iit.uncon_effect_repertoire(tpm)

ans = cei(subset, state, tpm);

print(ans);

#%%

