# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from __future__ import print_function
import json
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import numpy as np
from di2 import di2

def describe_sequences():
    aa_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    di_letters = ["%s%s" % (a, b) for a in aa_letters for b in aa_letters]
    letters = {1 : aa_letters, 2 : di_letters}
    
    def counter(string, seq_type):
        '''
        A function for counting the number of letters present.
        
        Returns a list of (letter, #occurances) tuples. 
        '''
        l = len(string)
        d = {i : 0 for i in letters[seq_type]}
        if seq_type == 1:
            for s in string:
                try:
                    d[s] += 1.0
                except KeyError:
                    d[s] = 1.0
            d = {k : d[k]/l for k in d}
        if seq_type == 2:        
            for a in range(l-1):
                s = string[a:a+seq_type]
                try:
                    d[s] += 1.0
                except KeyError:
                    d[s] = 1.0
            d = {k : d[k]/(l-1) for k in d}
        return d
        
    def residue_distribution(all_residues, seq_type):
        '''
        Takes as arguments a string with letters, and the type of sequence represented.
        Returns an alphabetically ordered string of relative frequencies, correct to three decimal places. 
        '''
        d = counter(all_residues, seq_type)
        residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] ]))                              ##Removes ambiguous letters
        r_c = [i[1] for i in residue_counts]
        dis = np.array([r_c,])
        return dis
    
    peptides = [{"seq" : "FLPILASLAAKFGPKLFCLVTKKC", "cTer" : None, "activity" : "YES"},
                {"seq" : "ILGPVISTIGGVLGGLLKNL", "cTer" : "Amidation", "activity" : "YES"},
                {"seq": "GIGGKILSGLKTALKGAAKELASTYLH", "cTer" : None, "activity" : "NO"},
                {"seq": "GIGSAILSAGKSALKGLAKGLAEHFAN", "cTer" : None, "activity" : "NO"},
                {"seq": "FLSLIPHAINAVSAIAKHF", "cTer" : "Amidation", "activity" : "NO"},
    ]
    
    
    for peptide in peptides:
        #print(peptide["id"])
        #print(peptide["seq"])
        
        globdesc = GlobalDescriptor(peptide["seq"])
        globdesc.calculate_all(amide = peptide["cTer"] == "Amidation")
        
        #peptide["GlobalDescriptor"] = globdesc
        
        #print(peptide["GlobalDescriptor"].descriptor)
        
        #Eisenberg hydrophobicity consensus
        #Take most of the values from here
        
        pepdesc = PeptideDescriptor(peptide["seq"], "eisenberg")
        pepdesc.calculate_global()
        pepdesc.calculate_moment(append=True)
        #pepdesc.calculate_profile(append=True, prof_type = "uH")
        
        pepdesc.load_scale("Ez")
        pepdesc.calculate_global(append=True)
        
        pepdesc.load_scale("charge_phys")
        pepdesc.calculate_moment(append=True)
        pepdesc.calculate_global(append=True)
        
        pepdesc.load_scale("flexibility")
        pepdesc.calculate_moment(append=True)
        pepdesc.calculate_global(append=True)
        
        pepdesc.load_scale("polarity")
        pepdesc.calculate_moment(append=True)
        pepdesc.calculate_global(append=True)
        
        pepdesc.load_scale("isaeci")
        pepdesc.calculate_global(append=True)
    
        pepdesc.load_scale("refractivity")
        pepdesc.calculate_moment(append=True)
        pepdesc.calculate_global(append=True)
        
        pepdesc.load_scale("z5")
        pepdesc.calculate_global(append=True)
        
        #peptide["PeptideDescriptor"] = pepdesc
    
        peptide["TotalDescriptor"] = str(np.concatenate((pepdesc.descriptor, globdesc.descriptor), axis=1))
        
        try:
            pepid = np.array([[int(peptide["id"].replace("HEMOLYTIK",""))]])
        except KeyError:
            pepid = np.array([[0]])
        
        freq_1d = residue_distribution(peptide["seq"], 1)
        freq_2d = residue_distribution(peptide["seq"], 2)
        
        len_peptide = np.array([[len(peptide["seq"])]])
        
        if peptide["activity"] == "YES":
            pepact = 1
        else:
            pepact = 0
        pepact = np.array([[pepact]])
        
        peptide_di2 = di2(peptide["seq"])
        
        peptide["array"] = np.concatenate((pepid, pepdesc.descriptor, globdesc.descriptor, len_peptide, 
               freq_1d, 
               #freq_2d, 
               #peptide_di2, 
               pepact,), axis=1)
        #print(peptide["TotalDescriptor"])
        
    
    x = np.concatenate([peptide["array"] for peptide in peptides], axis=0)
    print(x)
    
    np.save("hemolytik_array_custom_tests", x, allow_pickle=False)
    
describe_sequences()