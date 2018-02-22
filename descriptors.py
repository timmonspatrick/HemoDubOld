# -*- coding: utf-8 -*-

from __future__ import print_function
import json
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import numpy as np
from di2 import di2, di3
from pydpi.pypro import PyPro
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from aaf import aaf
from PyBioMed.PyProtein import CTD, ConjointTriad
import random
from sklearn.model_selection import train_test_split

random.seed(31)

print("beginning")
def describe_sequences():
    path = r"C:\Users\Patrick\OneDrive - University College Dublin\Bioinformatics\HemolyticStudies\BOTH_peptides.json"

    aa_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    di_letters = ["%s%s" % (a, b) for a in aa_letters for b in aa_letters]
    tri_letters = ["%s%s%s" % (a, b, c) for a in aa_letters for b in aa_letters for c in aa_letters]
    conjoint_letters = ["A", "I", "Y", "H", "R", "D", "C"]
    letters = {1 : aa_letters, 2 : di_letters, 3 : tri_letters, 4 : conjoint_letters}
    
    #Conjoint src = https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0828-1
    
    
    conjoint_dict = {"A" : "A", "G" : "A", "V" : "A",
                     "I" : "I", "L" : "I", "F" : "I", "P" : "I",
                     "Y" : "Y", "M" : "Y", "T" : "Y", "S" : "Y",
                     "H" : "H", "N" : "H", "Q" : "H", "W" : "H",
                     "R" : "R", "K" : "R",
                     "D" : "D", "E" : "D",
                     "C" : "C",}
    
    
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
        if seq_type == 3:        
            for a in range(l-2):
                s = string[a:a+seq_type]
                try:
                    d[s] += 1.0
                except KeyError:
                    d[s] = 1.0
            d = {k : d[k]/(l-2) for k in d}
        return d
    
    def counter_boolean(string, seq_type):
        '''
        A function for counting the number of letters present.
        
        Returns a list of (letter, #occurances) tuples. 
        '''
        l = len(string)
        d = {i : 0 for i in letters[seq_type]}
        if seq_type == 1:
            for s in string:
                try:
                    d[s] = 1.0
                except KeyError:
                    d[s] = 1.0
        if seq_type == 2:        
            for a in range(l-1):
                s = string[a:a+seq_type]
                try:
                    d[s] = 1.0
                except KeyError:
                    d[s] = 1.0
        return d
    
    def counter_abs(string, seq_type):
        '''
        A function for counting the number of letters present.
        
        Returns a list of (letter, #occurances) tuples. 
        '''
        l = len(string)
        d = {i : 0 for i in letters[seq_type]}
        if seq_type == 1:
            for s in string:
                try:
                    d[s] = d[s] + 1.0
                except KeyError:
                    d[s] = 1.0
        if seq_type == 2:        
            for a in range(l-1):
                s = string[a:a+seq_type]
                try:
                    d[s] = d[s] + 1.0
                except KeyError:
                    d[s] = 1.0
        return d
        
    def residue_distribution(all_residues, seq_type, dp):
        '''
        Takes as arguments a string with letters, and the type of sequence represented.
        Returns an alphabetically ordered string of relative frequencies, correct to three decimal places. 
        '''
        d = counter(all_residues, seq_type)
        if seq_type == 1:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] ]))                              ##Removes ambiguous letters
        elif seq_type == 2:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] if dp[i] >= 50])) 
        elif seq_type == 3:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] if tp[i] >= 20])) 
        elif seq_type == 4:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] ]))               

        r_c = [i[1] for i in residue_counts]
        dis = np.array([r_c,])
        return dis
    
    def residue_boolean(all_residues, seq_type, dp):
        '''
        Takes as arguments a string with letters, and the type of sequence represented.
        Returns an alphabetically ordered string of relative frequencies, correct to three decimal places. 
        '''
        d = counter_boolean(all_residues, seq_type)
        if seq_type == 1:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] ]))                              ##Removes ambiguous letters
        elif seq_type == 2:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] if dp[i] >= 50])) 
        r_c = [i[1] for i in residue_counts]
        dis = np.array([r_c,])
        return dis
    
    def residue_abs(all_residues, seq_type, dp):
        '''
        Takes as arguments a string with letters, and the type of sequence represented.
        Returns an alphabetically ordered string of relative frequencies, correct to three decimal places. 
        '''
        d = counter_abs(all_residues, seq_type)
        if seq_type == 1:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] ]))                              ##Removes ambiguous letters
        elif seq_type == 2:
            residue_counts = list(sorted([(i, d[i]) for i in letters[seq_type] if dp[i] >= 50])) 
        r_c = [i[1] for i in residue_counts]
        dis = np.array([r_c,])
        return dis
    
    with open(path, "r") as f:
        text = f.read()
    
    peptides = eval(text)["Peptides"]
    
    train_peptides, test_peptides = train_test_split(peptides, test_size=0.15, random_state=42)
    
    train_peptides_seqs = [peptide["seq"] for peptide in train_peptides]
    
    for peptide in peptides:
        if peptide["seq"] in train_peptides_seqs:
            peptide["train"] = True
        else:
            peptide["train"] = False
    
    print(len([p for p in peptides if p["train"] == True]))
    print(len([p for p in peptides if p["train"] == False]))
    
    new_peptides = []
    for peptide in peptides:
        if peptide["train"] == True:
            new_peptide = peptide.copy()
            new_seq = ''.join(reversed(peptide["seq"]))
            new_peptide["seq"] = new_seq
            new_peptides.append(new_peptide)
            
    #peptides.extend(new_peptides)
    random.shuffle(peptides)


    print(len([p for p in peptides if p["train"] == True]))
    print(len([p for p in peptides if p["train"] == False]))
    print("doubling complete")
    
    dp = {i: 0 for i in letters[2]}
    tp = {i: 0 for i in letters[3]}
    
    name_i = 0

    
    for peptide in peptides:
        temp_set = set()
        seq = peptide["seq"]
        l = len(seq)
        for a in range(l-1):
            s = seq[a:a+2]
            temp_set.add(s)
        for s in temp_set:
            dp[s] = dp[s] + 1
    
    for peptide in peptides:
        temp_set = set()
        seq = peptide["seq"]
        l = len(seq)
        for a in range(l-2):
            s = seq[a:a+3]
            temp_set.add(s)
        for s in temp_set:
            tp[s] = tp[s] + 1
            
    for peptide in peptides:
        peptide["conjoint_seq"] = "".join([conjoint_dict[letter] for letter in peptide["seq"]])
    
    
    for peptide in peptides:
        
        
        
        globdesc = GlobalDescriptor(peptide["seq"])
        globdesc.calculate_all(amide = peptide["cTer"] == "Amidation")
        
        ctdc = CTD.CalculateC(peptide["seq"])
        ctdc_keys = list(sorted(list([key for key in ctdc])))
        ctdc_vals = np.array([[ctdc[key] for key in ctdc_keys]]) 
        
        conjointtriad = ConjointTriad.CalculateConjointTriad(peptide["seq"])
        conjointtriad_keys = list(sorted(list([key for key in conjointtriad])))
        conjointtriad_vals = np.array([[conjointtriad[key] for key in conjointtriad_keys]]) 
        
        conjoint_dis = residue_distribution(peptide["conjoint_seq"], 4, None)

        #peptide["GlobalDescriptor"] = globdesc
        
        #print(peptide["GlobalDescriptor"].descriptor)
        
        #Eisenberg hydrophobicity consensus
        #Take most of the values from here
        
        pepdesc = PeptideDescriptor(peptide["seq"], "eisenberg")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        #pepdesc.calculate_profile(append=True, prof_type = "uH")
        
        pepdesc.load_scale("Ez")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("aasi")
        pepdesc.calculate_global(append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("abhprk")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("charge_acid")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("cougar")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("gravy")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("hopp-woods")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("kytedoolittle")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("ppcali")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("msw")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("charge_phys")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("flexibility")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("bulkiness")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("TM_tend")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("mss")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("t_scale")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("peparc")
        pepdesc.calculate_arc(modality="max", append=True)
        pepdesc.calculate_arc(modality="mean", append=True)        

        pepdesc.load_scale("msw")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("polarity")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("pepcats")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("isaeci")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
    
        pepdesc.load_scale("refractivity")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("z3")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("z5")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)

        #pepdesc.load_scale("PPCALI")
        #pepdesc.calculate_autocorr(2)        
        #peptide["PeptideDescriptor"] = pepdesc
        
        protein = PyPro()
        protein.ReadProteinSequence(peptide["seq"])
        paac=protein.GetPAAC(lamda=1, weight=0.05)
        paac2 = [[
                paac[a] for a in list(sorted([k for k in paac], key = lambda x : int(x.replace("PAAC",""))))
                ]]
        
        cTer = np.array([[1 if peptide["cTer"] == "Amidation" else 0]])
        
        
        paac = np.array(paac2)
        
        analysed_seq = ProteinAnalysis(peptide["seq"])
        secondary_structure_fraction = np.array([analysed_seq.secondary_structure_fraction()])
        
        
        peptide["TotalDescriptor"] = str(np.concatenate((pepdesc.descriptor, globdesc.descriptor), axis=1))
        
        try:
            pepid = np.array([[int(peptide["id"].replace("HEMOLYTIK","").replace("DRAMP","").replace("DBAASP",""))]])
        except KeyError:
            pepid = 0
            
        pep_train = np.array([[1 if peptide["train"] == True else 0]])
        
        freq_1d = residue_distribution(peptide["seq"], 1, dp)
        freq_2d = residue_distribution(peptide["seq"], 2, dp)
        freq_3d = residue_distribution(peptide["seq"], 3, dp)
        freq_1dbool = residue_boolean(peptide["seq"], 1, dp)
        freq_2dbool = residue_boolean(peptide["seq"], 2, dp)
        freq_1dabs = residue_abs(peptide["seq"], 1, dp)
        freq_2dabs = residue_abs(peptide["seq"], 2, dp)
        
        len_peptide = np.array([[len(peptide["seq"])]])
        
        if peptide["activity"] == "YES":
            pepact = 1
        else:
            pepact = 0
        pepact = np.array([[pepact]])
        
        peptide_di2 = di2(peptide["seq"])
        peptide_di3 = di3(peptide["conjoint_seq"])
        
        ####################### AAindex #########################
        to_get = [("CHAM810101", "mean"), #Steric Hinderance
                ("CHAM810101", "total"), #Steric Hinderance
                  ("KYTJ820101", "mean"), #Hydropathy
                  ("KLEP840101", "total"), #Charge
                  ("KLEP840101", "mean"), #Charge
                  ("MITS020101", "mean"), #Amphiphilicity
                  ("FAUJ830101", "mean"), #Hydrophobic parameter pi
                  ("GOLD730102", "total"), #Residue volume
                  ("MEEJ800101", "mean"), #Retention coefficient in HPLC
                  ("OOBM850105", "mean"), #Optimized side chain interaction parameter
                  ("OOBM850105", "total"), #Optimized side chain interaction parameter
                  ("VELV850101", "total"), #Electron-ion interaction parameter
                  ("VELV850101", "mean"), #Electron-ion interaction parameter
                  ("PUNT030102", "mean"), #Knowledge-based membrane-propensity scale from 3D_Helix
                  ("BHAR880101", "mean"), #Average flexibility indeces
                  ("KRIW790102", "mean"), #Fraction of site occupied by water
                  ("PLIV810101", "mean"), #Partition coefficient
                  ("ZIMJ680102", "mean"), #Bulkiness
                  ("ZIMJ680102", "total"), #Bulkiness
                  ("ZHOH040101", "mean"), #Stability scale
                  ("CHAM820102", "total"), #Free energy solubility in water
                                          #From HemoPi: src = https://github.com/riteshcanfly/Hemopi/blob/master/pcCalculator.java
                  ("HOPT810101", "mean"), #Hydrophilicity 
                  ("EISD840101", "mean"), #Hydrophobicity
                  ("FAUJ880109", "total"), #Net Hydrogen
                  ("EISD860101", "mean"), #Solvation
                ]
        
        tetra_peptides = ["KLLL",    # src = https://github.com/riteshcanfly/Hemopi/blob/master/tetrapos.txt
                          "GCSC",
                          "AAAK",
                          "KLLS",
                          "LGKL",
                          "VLKA",
                          "LLGK",
                          "LVGA",
                          "LSDF",
                          "SDFK",
                          "SWLR",
                          "WLRD",]
        
        tp_bin = [] 
        for t_p in tetra_peptides:
            if t_p in peptide["seq"]:
                tp_bin.append(1)
            else:
                tp_bin.append(0)
        tp_bin = np.array([tp_bin])
        
        for identifier, mode in to_get:
            x = aaf(peptide["seq"], identifier, mode)
        
        aminoacidindeces = np.array([[aaf(peptide["seq"], identifier, mode) for identifier, mode in to_get]])
            
        peptide["array"] = np.concatenate((pepid, pep_train, pepdesc.descriptor, globdesc.descriptor, len_peptide, 
               cTer,
               secondary_structure_fraction,
               aminoacidindeces,
               ctdc_vals,
               conjointtriad_vals, 
               tp_bin, 
               freq_1d, 
               freq_2d,
               freq_3d,
               freq_1dbool,
               freq_2dbool,
               freq_1dabs,
               freq_2dabs,
               peptide_di2,
               peptide_di3,  #Conjoint Alphabet
               paac,
               pepact,), axis=1)
        #print(peptide["TotalDescriptor"])
    
    
    x = np.concatenate([peptide["array"] for peptide in peptides], axis=0)
    
    np.save("peptides_array", x, allow_pickle=False)
    
describe_sequences()
print("complete")