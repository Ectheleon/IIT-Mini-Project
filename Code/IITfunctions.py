#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 15:46:28 2017

@author: Jonathan Peters
"""

import IIT2 as iit
import numpy as np

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
    
    print(f_uncon)
    print(p_uncon)
    
    f = iit.effect_repertoire(subset, full_set, state, tpm, base);
    p = iit.cause_repertoire(subset, full_set, state, tpm, base);
    
    print(f)
    print(p)
    
    cause_information = iit.EMD1(p_uncon, p);
    effect_information = iit.EMD1(f_uncon, f);
    
    print('ci',cause_information)
    print('ei',effect_information)
    
    return np.minimum(cause_information, effect_information);



#%%

subset = set([0]);
state = 1;
tpm = np.array([[0,0,0],[1,0,0],[1,0,1],[1,0,1],[0,0,1], [1,1,1], [1,0,0],[1,1,0]]);


iit.uncon_effect_repertoire(tpm)

ans = cei(subset, state, tpm);

print(ans);