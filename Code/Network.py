#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:08:37 2017

@author: Jonathan Peters
"""

#I'm not yet sure why I'll need the node class, but first I need to match the 
#matlab code before I can improve it

import numpy as np

#%%
class Nodes(object):
    
    def __init__(self, nnodes, tpm):
        self.nnodes = nnodes;
        self.num = 2*nnodes;
        
#%%


class Network(object):
    
    def __init__(self, nnodes, tpm,past, present, options, connectivity=False)):
        self.nnodes = nnodes;
        self.full_system = range(nnodes);
        self.nsubsets = 2**nnodes;
        
        if connectivity:
            self.connectivity_matrix = connectivity
        else:
            self.connectivity_matrix = np.ones([nnodes, nnodes]);
        
        self.tpm = tpm;
        self.past = past;
        self.present = present;
        self.options = options;
        
    
    
        
        
    