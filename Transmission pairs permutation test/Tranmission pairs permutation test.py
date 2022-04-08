#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 23:54:48 2021

@author: Chong Li
"""
import random
import collections
import pandas as pd
import seaborn as sns
from random import sample
from collections import Counter

with open("Nasal_H3.csv", "r") as Nasal_H3:
     
     
     Mutation_H3=[]
     Treatment_H3=[]
     PigID_H3=[]
     SampleID_H3=[]
     Freq_H3=[]
     title_H3=Nasal_H3.readline()
     
     for line_H3 in Nasal_H3:
         line_H3=line_H3.strip()
         Nasal_H3=line_H3.split(sep=',')
         Mutation_H3.append(Nasal_H3[2]+Nasal_H3[1]+Nasal_H3[18]+Nasal_H3[23])
         Treatment_H3.append(Nasal_H3[21])
         PigID_H3.append(Nasal_H3[19])
         SampleID_H3.append(Nasal_H3[0])
         Freq_H3.append(Nasal_H3[6])



H3_variants={}

H3_freq={}


for index, Pig in enumerate (PigID_H3):
    if Pig not in H3_variants.keys():
       H3_variants[Pig]={SampleID_H3[index]:[]}
       
    else:
       H3_variants[Pig][SampleID_H3[index]]=[] 

for No, Sample in enumerate (SampleID_H3): 
    H3_variants[PigID_H3[No]][Sample].append(Mutation_H3[No])
    
    
for index, Pig in enumerate (PigID_H3):
    if Pig not in H3_freq.keys():
       H3_freq[Pig]={SampleID_H3[index]:[]}
       
    else:
       H3_freq[Pig][SampleID_H3[index]]=[] 

for No, Sample in enumerate (SampleID_H3): 
    H3_freq[PigID_H3[No]][Sample].append(float(Freq_H3[No]))
    
    
def shared_variant_proportion (pig_1, sample_1, pig_2, sample_2, variant_dict):
    
    variant_pool=variant_dict[pig_1][sample_1]+variant_dict[pig_2][sample_2]
            
    count=collections.Counter(variant_pool).values()
            
    shared_number = 0
    for i in count:
        if int(i) > 1:
           shared_number += 1
    total = int(len(variant_dict[pig_1][sample_1]))+ int(len(variant_dict[pig_2][sample_2]))      
    proportion = int(shared_number)*2/ total
    return proportion  

#def pumutation_random_pairs (variant_dict, stimulation):
    stimulated_proportion=[]
    for a in range(0, stimulation, 1):
        sample_pair_x = random.sample(variant_dict.keys(), 1)
        pig_x = sample_pair_x[0]
        sample_pair_y = random.sample(variant_dict.keys(), 1)
        pig_y = sample_pair_y[0]
        sample_x = random.sample(variant_dict[pig_x].keys(),1)[0]
        sample_y = random.sample(variant_dict[pig_y].keys(),1)[0]
        
        z=shared_variant_proportion (pig_x, sample_x, pig_y, sample_y, variant_dict)
        
        stimulated_proportion.append(float(z))
        
    return stimulated_proportion

def pumutation_random_pairs (variant_dict, stimulation):
    stimulated_proportion=[]
    for a in range(0, stimulation, 1):
        sample_pair = random.sample(variant_dict.keys(), 2)
        pig_x = sample_pair[0]
        
        pig_y = sample_pair[1]
        sample_x = random.sample(variant_dict[pig_x].keys(),1)[0]
        sample_y = random.sample(variant_dict[pig_y].keys(),1)[0]
        z=shared_variant_proportion (pig_x, sample_x, pig_y, sample_y, variant_dict)
        
        stimulated_proportion.append(float(z))
        
    return stimulated_proportion


# run the test
stimulation = 100000
H3_nasal = pumutation_random_pairs (H3_variants, stimulation)


df_H3_nasal = pd.DataFrame({"Shared proportion": H3_nasal})
df_H3_nasal = df_H3_nasal.reset_index()
df_H3_nasal["Group"]= "Random pair" 



df_transmission_permutation=pd.concat([df_H3_nasal])
           
df_transmission_permutation.to_csv(r'/Users/apple/Desktop/H3_Permutation_pairs.csv',index=False)     
    