# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 15:02:19 2017

@author: Patrick
"""
from __future__ import print_function
import numpy as np
conjoint_letters = ["A", "I", "Y", "H", "R", "D", "C"]
aa_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
di_letters = ["%s%s" % (a, b) for a in aa_letters for b in aa_letters]
di3_letters = ["%s%s%s" % (a, b, c) for a in aa_letters for b in aa_letters for c in conjoint_letters]

def counter(string):
        '''
        A function for counting the number of letters present.
        
        Returns a list of (letter, #occurances) tuples. 
        '''
        
        l = max(1, len(string))
        d = {i : 0 for i in di_letters}       
        for s in string:
            try:
                d[s] += 1.0
            except KeyError:
                d[s] = 1.0
        d = {k : d[k]/(l) for k in d}
        return d
    
def counter3(string):
        '''
        A function for counting the number of letters present.
        
        Returns a list of (letter, #occurances) tuples. 
        '''
        
        l = max(1, len(string))
        d = {i : 0 for i in conjoint_letters}       
        for s in string:
            try:
                d[s] += 1.0
            except KeyError:
                d[s] = 1.0
        d = {k : d[k]/(l) for k in d}
        return d
        
def residue_distribution(all_residues):
    '''
    Takes as arguments a string with letters, and the type of sequence represented.
    Returns an alphabetically ordered string of relative frequencies, correct to three decimal places. 
    '''
    d = counter(all_residues)
    residue_counts = list(sorted([(i, d[i]) for i in di_letters ]))                              ##Removes ambiguous letters
    r_c = [i[1] for i in residue_counts]
    dis = np.array([r_c,])
    return dis

def residue_distribution3(all_residues):
    '''
    Takes as arguments a string with letters, and the type of sequence represented.
    Returns an alphabetically ordered string of relative frequencies, correct to three decimal places. 
    '''
    d = counter3(all_residues)
    residue_counts = list(sorted([(i, d[i]) for i in conjoint_letters ]))                              ##Removes ambiguous letters
    r_c = [i[1] for i in residue_counts]
    dis = np.array([r_c,])
    return dis

def di2(seq):
    '''
    A function to return all the di2s for a sequence.
    Eg. ABCDEF --> AD, BE, CF
    '''
    l = []
    for a in range(len(seq)):
        try:
            x = "%s%s" % (seq[a], seq[a + 3 ])
            l.append(x)
        except IndexError:
            pass
    
    return residue_distribution(l)

def di3(seq):
    '''
    A function to return all the di3s for a sequence.
    Eg. ABCDEFGHI --> ADG, BEH, CFI
    '''
    l = []
    for a in range(len(seq)):
        try:
            x = "%s%s%s" % (seq[a], seq[a + 3 ], seq[a + 6])
            l.append(x)
        except IndexError:
            pass
    
    return residue_distribution(l)
    
    