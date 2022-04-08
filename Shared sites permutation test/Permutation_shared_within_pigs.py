#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 18:18:50 2021

@author: Chong Li
"""
import random
import collections
import pandas as pd
import seaborn as sns
from random import sample
from collections import Counter

    
#Generate the snp file for H1 virus.    

with open("SNP_H1_protein.csv", "r") as SNP_H1_protein:
    

    Sample_H1=[]
    Group_H1=[]
    Room_H1=[]
    Pig_H1=[]
    PB2_H1=[]
    PB1_H1=[]
    PB1_F2_H1=[]
    PA_H1=[]
    PA_X_H1=[]
    NP_H1=[]
    HA_H1=[]
    NA_H1=[]
    M1_H1=[]
    M2_H1=[]
    NS1_H1=[]
    NS2_H1=[]
    title_H1=SNP_H1_protein.readline()

    for line_H1 in SNP_H1_protein:
        line_H1=line_H1.strip()
        split_file_H1=line_H1.split(sep=',')
        Sample_H1.append(split_file_H1[0])
        Group_H1.append(split_file_H1[2])
        Room_H1.append(split_file_H1[4])
        Pig_H1.append(split_file_H1[5])
        PB2_H1.append(split_file_H1[6])
        PB1_H1.append(split_file_H1[7])
        PB1_F2_H1.append(split_file_H1[8])
        PA_H1.append(split_file_H1[9])
        PA_X_H1.append(split_file_H1[10])
        HA_H1.append(split_file_H1[11])
        NP_H1.append(split_file_H1[12])
        NA_H1.append(split_file_H1[13])
        M1_H1.append(split_file_H1[14])
        M2_H1.append(split_file_H1[15])
        NS1_H1.append(split_file_H1[16])
        NS2_H1.append(split_file_H1[17])

PRIME_BOOST_H1={}
SINGLE_LAIV_H1={}
NO_VAC_H1={}

for index_H1, words_H1 in enumerate (Group_H1):
    if words_H1 == "PRIME BOOST":
       if Pig_H1[int(index_H1)] not in PRIME_BOOST_H1.keys():       
          PRIME_BOOST_H1[Pig_H1[int(index_H1)]] = {Sample_H1[int(index_H1)]:{"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}} 
       else:
          PRIME_BOOST_H1[Pig_H1[int(index_H1)]][Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]} 
       
    if words_H1 == "SINGLE LAIV":
       if Pig_H1[int(index_H1)] not in SINGLE_LAIV_H1.keys():
          SINGLE_LAIV_H1[Pig_H1[int(index_H1)]] = {Sample_H1[int(index_H1)]:{"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}} 
       else:
          SINGLE_LAIV_H1[Pig_H1[int(index_H1)]][Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}
    
    if words_H1 == "NO VAC":
       if Pig_H1[int(index_H1)] not in NO_VAC_H1.keys():
          NO_VAC_H1[Pig_H1[int(index_H1)]] = {Sample_H1[int(index_H1)]:{"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}} 
       else:
          NO_VAC_H1[Pig_H1[int(index_H1)]][Sample_H1[int(index_H1)]]= {"PB2":[PB2_H1[int(index_H1)]], "PB1":[PB1_H1[int(index_H1)]], "PB1-F2":[PB1_F2_H1[int(index_H1)]], "PA":[PA_H1[int(index_H1)]], "PA-X":[PA_X_H1[int(index_H1)]], "HA":[HA_H1[int(index_H1)]], "NP":[NP_H1[int(index_H1)]],"NA":[NA_H1[int(index_H1)]], "M1":[M1_H1[int(index_H1)]], "M2":[M2_H1[int(index_H1)]], "NS1":[NS1_H1[int(index_H1)]], "NS2":[NS2_H1[int(index_H1)]]}

# Generate the snp files for H3 virus.
with open("SNP_H3_protein.csv", "r") as SNP_H3_protein:
    

    Sample_H3=[]
    Group_H3=[]
    Room_H3=[]
    Pig_H3=[]
    PB2_H3=[]
    PB1_H3=[]
    PB1_F2_H3=[]
    PA_H3=[]
    PA_X_H3=[]
    NP_H3=[]
    HA_H3=[]
    NA_H3=[]
    M1_H3=[]
    M2_H3=[]
    NS1_H3=[]
    NS2_H3=[]
    title_H3=SNP_H3_protein.readline()

    for line_H3 in SNP_H3_protein:
        line_H3=line_H3.strip()
        split_file_H3=line_H3.split(sep=',')
        Sample_H3.append(split_file_H3[0])
        Group_H3.append(split_file_H3[2])
        Room_H3.append(split_file_H3[4])
        Pig_H3.append(split_file_H3[5])
        PB2_H3.append(split_file_H3[6])
        PB1_H3.append(split_file_H3[7])
        PB1_F2_H3.append(split_file_H3[8])
        PA_H3.append(split_file_H3[9])
        PA_X_H3.append(split_file_H3[10])
        HA_H3.append(split_file_H3[11])
        NP_H3.append(split_file_H3[12])
        NA_H3.append(split_file_H3[13])
        M1_H3.append(split_file_H3[14])
        M2_H3.append(split_file_H3[15])
        NS1_H3.append(split_file_H3[16])
        NS2_H3.append(split_file_H3[17])

PRIME_BOOST_H3={}
SINGLE_LAIV_H3={}
NO_VAC_H3={}

for index_H3, words_H3 in enumerate (Group_H3):
    if words_H3 == "PRIME BOOST":
        
       if Pig_H3[int(index_H3)] not in PRIME_BOOST_H3.keys(): 
          PRIME_BOOST_H3[Pig_H3[int(index_H3)]]= {Sample_H3[int(index_H3)]:{"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]}}
       else:
          PRIME_BOOST_H3[Pig_H3[int(index_H3)]][Sample_H3[int(index_H3)]] = {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 
    
    if words_H3 == "SINGLE LAIV":
       if Pig_H3[int(index_H3)] not in SINGLE_LAIV_H3.keys(): 
          SINGLE_LAIV_H3[Pig_H3[int(index_H3)]]= {Sample_H3[int(index_H3)]:{"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]}}
       else:
          SINGLE_LAIV_H3[Pig_H3[int(index_H3)]][Sample_H3[int(index_H3)]] = {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 

    if words_H3 == "NO VAC":
       if Pig_H3[int(index_H3)] not in NO_VAC_H3.keys(): 
          NO_VAC_H3[Pig_H3[int(index_H3)]]= {Sample_H3[int(index_H3)]:{"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]}}
       else:
          NO_VAC_H3[Pig_H3[int(index_H3)]][Sample_H3[int(index_H3)]] = {"PB2":[PB2_H3[int(index_H3)]], "PB1":[PB1_H3[int(index_H3)]], "PB1-F2":[PB1_F2_H3[int(index_H3)]], "PA":[PA_H3[int(index_H3)]], "PA-X":[PA_X_H3[int(index_H3)]], "HA":[HA_H3[int(index_H3)]], "NP":[NP_H3[int(index_H3)]],"NA":[NA_H3[int(index_H3)]], "M1":[M1_H3[int(index_H3)]], "M2":[M2_H3[int(index_H3)]], "NS1":[NS1_H3[int(index_H3)]], "NS2":[NS2_H3[int(index_H3)]]} 

#Generate the set of the numbers of nucleotides in each gene segment.
Nucleotide_region_H1={}
Nucleotide_region_H3={}

for num_H1 in range(0,len(Pig_H1),1):
    Nucleotide_region_H1[Pig_H1[int(num_H1)]]={"PB2":759, "PB1":757, "PB1-F2":79, "PA":716, "PA-X":232, "HA":566, "NP":498, "NA":469, "M1":252, "M2":97, "NS1":219, "NS2":121}
for num_H3 in range(0,len(Pig_H3),1):
    Nucleotide_region_H3[Pig_H3[int(num_H3)]]={"PB2":759, "PB1":757, "PB1-F2":79, "PA":716, "PA-X":232, "HA":566, "NP":498, "NA":469, "M1":252, "M2":97, "NS1":219, "NS2":121}

#Generate the set of the 70% numbers of nucleotides in each gene segment.
Nucleotide_region_70_H1={}
Nucleotide_region_70_H3={}

for num_H1_70 in range(0,len(Pig_H1),1):
    Nucleotide_region_70_H1[Pig_H1[int(num_H1_70)]]={"PB2":531, "PB1":530, "PB1-F2":55, "PA":501, "PA-X":162, "HA":396, "NP":349, "NA":328, "M1":176, "M2":68, "NS1":153, "NS2":85}
for num_H3_70 in range(0,len(Pig_H3),1):
    Nucleotide_region_70_H3[Pig_H3[int(num_H3_70)]]={"PB2":531, "PB1":530, "PB1-F2":55, "PA":501, "PA-X":162, "HA":396, "NP":349, "NA":328, "M1":176, "M2":68, "NS1":153, "NS2":85}

#Generate the dictionary for pigID and room
Room_Pig_H1={"1":[], "3":[], "5":[], "7":[], "9":[], "11":[], "8":[], "10":[]}
Room_Pig_H3={"1":[], "3":[], "5":[], "7":[], "9":[], "11":[], "8":[], "10":[]}
for x in range(0,len(Room_H1),1):
    if Pig_H1[x] not in Room_Pig_H1[Room_H1[x]]:
       Room_Pig_H1[Room_H1[x]].append(Pig_H1[x]) 
       
for y in range(0,len(Room_H3),1):
    if Pig_H3[y] not in Room_Pig_H3[Room_H3[y]]:
       Room_Pig_H3[Room_H3[y]].append(Pig_H3[y])        
    
 
def permutation_test(SNP_array, Nucleotide_region, Room_dict, stimulation):
    def get_room_number(PigID):
        for Room in Room_dict:
            if PigID in Room_dict[Room]:
               return Room
 
        return "Not found."
    
    total_shared_sites=[]
    iteration_results = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[], "NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
    for a in range(0, stimulation, 1):
        temporary_results = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[],"NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
        temporary_results_test = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[],"NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
        for pigs in SNP_array:
            temporary_results_pigs = {"PB2":[], "PB1":[], "PB1-F2":[], "PA":[], "PA-X":[], "HA":[], "NP":[],"NA":[], "M1":[], "M2":[], "NS1":[], "NS2":[]}
            for samples in SNP_array[pigs]:
                for gene in SNP_array[pigs][samples]:
                    Nucleotides_range = range(0,Nucleotide_region[pigs][gene])
                    
                
                    Num_draw_list = SNP_array[pigs][samples][gene]
                    Num_draw = int(Num_draw_list[0])
                
            
                    random_draw = random.sample(Nucleotides_range, Num_draw)
                    
                    for r in random_draw:
                        if r not in temporary_results_pigs[gene]:
                           temporary_results_pigs[gene].append(r)
                           
            for gene in temporary_results_pigs:
                for h in temporary_results_pigs[gene]:
                    marker= str(h) +"-"+ get_room_number(pigs)
                    if marker not in temporary_results_test[gene]:
                       temporary_results_test[gene].append(marker) 
                       temporary_results[gene].append(h)
                    
        for gene in temporary_results:
            count=collections.Counter(temporary_results[gene]).values()
            
            more_than_once = 0
            for i in count:
                if i > 1:
                   more_than_once += 1
            
            iteration_results[gene].append(more_than_once)

    for c in range (0, stimulation, 1):
        shared_sites = iteration_results["PB2"][c]+iteration_results["PB1"][c]+ iteration_results["PB1-F2"][c]+iteration_results["PA"][c]+iteration_results["PA-X"][c]+iteration_results["HA"][c]+iteration_results["NP"][c]+iteration_results["NA"][c]+iteration_results["M1"][c]+iteration_results["M2"][c]+iteration_results["NS1"][c]+iteration_results["NS2"][c]
        total_shared_sites.append(shared_sites)

    return(total_shared_sites) 

#run the method.
stimulation = 100000
PRIMEBOOST_H1= permutation_test(PRIME_BOOST_H1, Nucleotide_region_H1, Room_Pig_H1, stimulation)
SINGLELAIV_H1= permutation_test(SINGLE_LAIV_H1, Nucleotide_region_H1, Room_Pig_H1, stimulation)
NOVAC_H1= permutation_test(NO_VAC_H1, Nucleotide_region_H1, Room_Pig_H1, stimulation)

PRIMEBOOST_H3= permutation_test(PRIME_BOOST_H3, Nucleotide_region_H3, Room_Pig_H3, stimulation)
SINGLELAIV_H3= permutation_test(SINGLE_LAIV_H3, Nucleotide_region_H3, Room_Pig_H3, stimulation)
NOVAC_H3= permutation_test(NO_VAC_H3, Nucleotide_region_H3, Room_Pig_H3, stimulation)
# run the method for 70% amino acid sites.
stimulation = 100000
PRIMEBOOST_H1_70= permutation_test(PRIME_BOOST_H1, Nucleotide_region_70_H1, Room_Pig_H1, stimulation)
SINGLELAIV_H1_70= permutation_test(SINGLE_LAIV_H1, Nucleotide_region_70_H1, Room_Pig_H1, stimulation)
NOVAC_H1_70= permutation_test(NO_VAC_H1, Nucleotide_region_70_H1, Room_Pig_H1, stimulation)

PRIMEBOOST_H3_70= permutation_test(PRIME_BOOST_H3, Nucleotide_region_70_H3, Room_Pig_H3, stimulation)
SINGLELAIV_H3_70= permutation_test(SINGLE_LAIV_H3, Nucleotide_region_70_H3, Room_Pig_H3, stimulation)
NOVAC_H3_70= permutation_test(NO_VAC_H3, Nucleotide_region_70_H3, Room_Pig_H3, stimulation)
# H1 convert dataset.
df_H1_PRIME_BOOST = pd.DataFrame({"whole_genome":PRIMEBOOST_H1})
df_H1_PRIME_BOOST = df_H1_PRIME_BOOST.reset_index()
df_H1_PRIME_BOOST ["Group"]= "PRIME BOOST"  

df_H1_SINGLE_LAIV = pd.DataFrame({"whole_genome":SINGLELAIV_H1})
df_H1_SINGLE_LAIV = df_H1_SINGLE_LAIV.reset_index()
df_H1_SINGLE_LAIV ["Group"]= "SINGLE LAIV"  

df_H1_NO_VAC = pd.DataFrame({"whole_genome":NOVAC_H1})
df_H1_NO_VAC = df_H1_NO_VAC.reset_index()
df_H1_NO_VAC ["Group"]= "NO VAC"  

df_H1_protein=pd.concat([df_H1_PRIME_BOOST,df_H1_SINGLE_LAIV,df_H1_NO_VAC])
#H3 convert dataset.
df_H3_PRIME_BOOST = pd.DataFrame({"whole_genome":PRIMEBOOST_H3})
df_H3_PRIME_BOOST = df_H3_PRIME_BOOST.reset_index()
df_H3_PRIME_BOOST ["Group"]= "PRIME BOOST"  

df_H3_SINGLE_LAIV = pd.DataFrame({"whole_genome":SINGLELAIV_H3})
df_H3_SINGLE_LAIV = df_H3_SINGLE_LAIV.reset_index()
df_H3_SINGLE_LAIV ["Group"]= "SINGLE LAIV"  

df_H3_NO_VAC = pd.DataFrame({"whole_genome":NOVAC_H3})
df_H3_NO_VAC = df_H3_NO_VAC.reset_index()
df_H3_NO_VAC ["Group"]= "NO VAC"  

df_H3_protein=pd.concat([df_H3_PRIME_BOOST,df_H3_SINGLE_LAIV,df_H3_NO_VAC])

# H1 convert dataset for 70% sites.
df_H1_PRIME_BOOST_70 = pd.DataFrame({"whole_genome":PRIMEBOOST_H1_70})
df_H1_PRIME_BOOST_70 = df_H1_PRIME_BOOST_70.reset_index()
df_H1_PRIME_BOOST_70 ["Group"]= "PRIME BOOST"  

df_H1_SINGLE_LAIV_70 = pd.DataFrame({"whole_genome":SINGLELAIV_H1_70})
df_H1_SINGLE_LAIV_70 = df_H1_SINGLE_LAIV_70.reset_index()
df_H1_SINGLE_LAIV_70 ["Group"]= "SINGLE LAIV"  

df_H1_NO_VAC_70 = pd.DataFrame({"whole_genome":NOVAC_H1_70})
df_H1_NO_VAC_70 = df_H1_NO_VAC_70.reset_index()
df_H1_NO_VAC_70 ["Group"]= "NO VAC"  

df_H1_70_protein=pd.concat([df_H1_PRIME_BOOST_70,df_H1_SINGLE_LAIV_70,df_H1_NO_VAC_70])

#H3 convert dataset for 70% sites.
df_H3_PRIME_BOOST_70 = pd.DataFrame({"whole_genome":PRIMEBOOST_H3_70})
df_H3_PRIME_BOOST_70 = df_H3_PRIME_BOOST_70.reset_index()
df_H3_PRIME_BOOST_70 ["Group"]= "PRIME BOOST"  

df_H3_SINGLE_LAIV_70 = pd.DataFrame({"whole_genome":SINGLELAIV_H3_70})
df_H3_SINGLE_LAIV_70 = df_H3_SINGLE_LAIV_70.reset_index()
df_H3_SINGLE_LAIV_70 ["Group"]= "SINGLE LAIV"  

df_H3_NO_VAC_70 = pd.DataFrame({"whole_genome":NOVAC_H3_70})
df_H3_NO_VAC_70 = df_H3_NO_VAC_70.reset_index()
df_H3_NO_VAC_70 ["Group"]= "NO VAC"  

df_H3_70_protein=pd.concat([df_H3_PRIME_BOOST_70,df_H3_SINGLE_LAIV_70,df_H3_NO_VAC_70])

#Export data.
df_H3_protein.to_csv(r'/Users/apple/Desktop/H3_Permutation_protein.csv',index=False)
df_H1_protein.to_csv(r'/Users/apple/Desktop/H1_Permutation_protein.csv',index=False)

#Export data for 70% sites.
df_H3_70_protein.to_csv(r'/Users/apple/Desktop/H3_Permutation_70%_protein.csv',index=False)
df_H1_70_protein.to_csv(r'/Users/apple/Desktop/H1_Permutation_70%_protein.csv',index=False)
