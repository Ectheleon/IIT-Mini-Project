# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

def buildTransition(rule_func, size):
    nstates = 2**size
    mat = np.array([rule_func(i) for i in range(nstates)])
    
    return mat        
        
        
def ex1(past_state):
    #This rule applies to the 3 node mechanism in the IIT paper. Node 1 is the 
    #'or' gate, 2 is the 'and' gate, and 3 is the 'xor' gate. I represent a 
    #state by an integer, where its binary form shows whether each gate is on
    #or off. i.e. 3 = 011 means node 1 is off, but the other two are on.
    
    if past_state > 7 or past_state < 0:
        raise ValueError('%n is not a valid state in a 3 node system', past_state)
    
    ans = 0
    
    if past_state & 3 >0:
        ans+=4
    if past_state &6>0 and past_state&6 <6:
        ans+=1
    if past_state & 5 ==5:
        ans+=2
    
    return ans

def repeat_rule(rule, n):
    #returns the transition array for the rule when applied n times. Requires 
    #the rule matrix as input
    nstates = np.size(rule);
    
    if n>1:
        new = np.array([rule[rule[i]] for i in range(nstates)])
        
        for i in range(0,n-2):
            print 'activated'
            new = np.array([new[rule[i]] for i in range(nstates)])
        
        return new;
    else:
        return rule
    
def unconFuture(rule):
    #finds the number of distinct states in the rule matrix, and calculates
    #the frequency of each distinct state
    nstates = np.size(rule)
    freqs = np.zeros(nstates)
    
    for i in range(nstates):
        freqs[rule[i]]+=1
    
    return np.divide(freqs, 1.0*nstates)

def Hamming(a,b):
    #computes a^b, and then counts how many nonzero bits are found in a^b
    tmp = a^b;
    count = 0
    while tmp > 0:
        if tmp&1:
            count+=1
        tmp >>=1;
    
    return count

def D(dist1, dist2):
    #it is assumed that dist1 and dist2 are vectors of the same length
    d = np.array([Hamming(i,j) for i in range(np.size(dist1)) for j in range(np.size(dist2))])
    return d