#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 18:31:00 2021

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
    

with open("Nasal_H1.csv", "r") as Nasal_H1:
     
     
     Mutation_H1=[]
     Treatment_H1=[]
     PigID_H1=[]
     SampleID_H1=[]
     Freq_H1=[]
     title_H1=Nasal_H1.readline()
     
     for line_H1 in Nasal_H1:
         line_H1=line_H1.strip()
         Nasal_H1=line_H1.split(sep=',')
         Mutation_H1.append(Nasal_H1[2]+Nasal_H1[1]+Nasal_H1[18]+Nasal_H1[23])
         Treatment_H1.append(Nasal_H1[21])
         PigID_H1.append(Nasal_H1[19])
         SampleID_H1.append(Nasal_H1[0])
         Freq_H1.append(Nasal_H1[6])



H1_variants={}

H1_freq={}


for index, Pig in enumerate (PigID_H1):
    if Pig not in H1_variants.keys():
       H1_variants[Pig]={SampleID_H1[index]:[]}
       
    else:
       H1_variants[Pig][SampleID_H1[index]]=[] 

for No, Sample in enumerate (SampleID_H1): 
    H1_variants[PigID_H1[No]][Sample].append(Mutation_H1[No])
    
    
for index, Pig in enumerate (PigID_H1):
    if Pig not in H1_freq.keys():
       H1_freq[Pig]={SampleID_H1[index]:[]}
       
    else:
       H1_freq[Pig][SampleID_H1[index]]=[] 

for No, Sample in enumerate (SampleID_H1): 
    H1_freq[PigID_H1[No]][Sample].append(float(Freq_H1[No]))



    
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


def sample_pairs (pig_n, pig_m, variant_dict):
    pig_n_m = []
    for a in variant_dict[pig_n].keys():
        for b in variant_dict[pig_m].keys():
            
            f = float(shared_variant_proportion(pig_n, a, pig_m, b, variant_dict)) 
                
            pig_n_m.append(f) 
            
    return pig_n_m    




def within_households (variant_dict, freq_dict, transmission_pairs):
    within_house=[]
    mark=[]
    for pig_w in variant_dict.keys():
        for pig_y in variant_dict.keys():
            if pig_w != pig_y and pig_w +"-"+pig_y not in mark:
               mark.append(pig_w +"-"+pig_y)
                
               if pig_w +"-"+pig_y in transmission_pairs:
                  
                  within_house += sample_pairs (pig_w, pig_y, variant_dict)
            elif pig_w != pig_y and pig_y +"-"+pig_w not in mark:
                 mark.append(pig_y +"-"+pig_w)         
    
                 if pig_y +"-"+pig_w in transmission_pairs:
                  
                    within_house += sample_pairs (pig_w, pig_y, variant_dict)          
   
    return within_house

def without_households (variant_dict, freq_dict, transmission_pairs):
    without_house=[]
    mark=[]
    for pig_w in variant_dict.keys():
        for pig_y in variant_dict.keys():
            if pig_w != pig_y and pig_w +"-"+pig_y not in mark:
               
                
               if pig_w +"-"+pig_y not in transmission_pairs:
                  
                  without_house += sample_pairs (pig_w, pig_y, variant_dict)
                
                  mark.append(pig_w +"-"+pig_y)
                               
                  mark.append(pig_y+"-"+pig_w)
    return without_house

def within_hosts (variant_dict, freq_dict):
    within_host=[]
    
    for pig_w in variant_dict.keys():
        mark=[]
        for sample_w in variant_dict[pig_w].keys():
            for sample_y in variant_dict[pig_w].keys():
                if sample_w != sample_y and sample_w + "-" + sample_y not in mark:
                   f = shared_variant_proportion(pig_w, sample_w, pig_w, sample_y, variant_dict) 
                   within_host.append(float(f))
                   mark.append(sample_w +"-"+sample_y)
                   mark.append(sample_y+"-"+sample_w)

    return within_host


H3_transmission_pairs_3 = []
H3_transmission_pairs_4 = ["5174-5179", "4490-4551", "4490-4931", "4490-5167", "4551-4931", "4551-5167", "4931-5167", "4479-4484", "4479-4933", "4479-4945", "4479-5166", "4484-4933", "4484-4945", "4484-5166", "4933-4945", "4933-5166", "4945-5166" ]
H1_transmission_pairs_3 = []
H1_transmission_pairs_4 = ["4490-4495", "4473-4483", "4481-4484"]

H3_within_house = within_households (H3_variants, H3_freq, H3_transmission_pairs_4)
H3_without_house = without_households (H3_variants, H3_freq, H3_transmission_pairs_3)
H3_within_hosts = within_hosts (H3_variants, H3_freq)

df_H3_within_house = pd.DataFrame({"Shared_proportion": H3_within_house})
df_H3_within_house = df_H3_within_house.reset_index()
df_H3_within_house["Group"]= "Household pairs" 

df_H3_without_house = pd.DataFrame({"Shared_proportion": H3_without_house})
df_H3_without_house = df_H3_without_house.reset_index()
df_H3_without_house["Group"]= "Random pairs" 

df_H3_within_hosts = pd.DataFrame({"Shared_proportion": H3_within_hosts})
df_H3_within_hosts = df_H3_within_hosts.reset_index()
df_H3_within_hosts["Group"]= "Individual pairs" 



H1_within_house = within_households (H1_variants, H1_freq, H1_transmission_pairs_4)
H1_without_house = without_households (H1_variants, H1_freq, H1_transmission_pairs_3)
H1_within_hosts = within_hosts (H1_variants, H1_freq)

df_H1_within_house = pd.DataFrame({"Shared_proportion": H1_within_house})
df_H1_within_house = df_H1_within_house.reset_index()
df_H1_within_house["Group"]= "Household pairs" 

df_H1_without_house = pd.DataFrame({"Shared_proportion": H1_without_house})
df_H1_without_house = df_H1_without_house.reset_index()
df_H1_without_house["Group"]= "Random pairs" 

df_H1_within_hosts = pd.DataFrame({"Shared_proportion": H1_within_hosts})
df_H1_within_hosts = df_H1_within_hosts.reset_index()
df_H1_within_hosts["Group"]= "Individual pairs" 

df_H3_within_out=pd.concat([df_H3_within_house,df_H3_without_house, df_H3_within_hosts])
df_H1_within_out=pd.concat([df_H1_within_house,df_H1_without_house, df_H1_within_hosts])           

df_H3_within_out.to_csv(r'/Users/apple/Desktop/H3 shared proportion pairs.csv',index=False)                  
df_H1_within_out.to_csv(r'/Users/apple/Desktop/H1 shared proportion pairs.csv',index=False)       





