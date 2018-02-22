# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 18:07:22 2017

@author: Patrick
"""

import aaindex
aaindex.init(path='.')

def aaf(sequence, identifier, mode):
    x = aaindex.get(identifier)
    
    total = 0
    total_list = []
    for aa in sequence:
        total = total + x.get(aa) 
        total_list.append(x.get(aa))
        
    if mode == "mean":
        return total / len(sequence)
    elif mode == "total":
        return total
    elif mode == "max":
        return max(total_list)
    
