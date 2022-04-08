#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 00:12:20 2021

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

H3_PRIME_BOOST={}
H3_SINGLE_LAIV={}
H3_NO_VAC={}
for a, b in enumerate (Treatment_H3):
    if b == "PRIME BOOST":
       if PigID_H3[a] not in H3_PRIME_BOOST.keys():
          H3_PRIME_BOOST[PigID_H3[a]]=[str(SampleID_H3[a])]
       
       else:
          if str(SampleID_H3[a]) not in H3_PRIME_BOOST[PigID_H3[a]]:
             H3_PRIME_BOOST[PigID_H3[a]].append(str(SampleID_H3[a]))  


    elif b == "SINGLE LAIV":
       if PigID_H3[a] not in H3_SINGLE_LAIV.keys():
          H3_SINGLE_LAIV[PigID_H3[a]]=[str(SampleID_H3[a])]
       
       else:
          if str(SampleID_H3[a]) not in H3_SINGLE_LAIV[PigID_H3[a]]:
             H3_SINGLE_LAIV[PigID_H3[a]].append(str(SampleID_H3[a])) 
             
    elif b == "NO VAC":
        if PigID_H3[a] not in H3_NO_VAC.keys():
          H3_NO_VAC[PigID_H3[a]]=[str(SampleID_H3[a])]
       
        else:
          if str(SampleID_H3[a]) not in H3_NO_VAC[PigID_H3[a]]:
             H3_NO_VAC[PigID_H3[a]].append(str(SampleID_H3[a]))




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

H1_PRIME_BOOST={}
H1_SINGLE_LAIV={}
H1_NO_VAC={}
for a, b in enumerate (Treatment_H1):
    if b == "PRIME BOOST":
       if PigID_H1[a] not in H1_PRIME_BOOST.keys():
          H1_PRIME_BOOST[PigID_H1[a]]=[str(SampleID_H1[a])]
       
       else:
          if str(SampleID_H1[a]) not in H1_PRIME_BOOST[PigID_H1[a]]:
             H1_PRIME_BOOST[PigID_H1[a]].append(str(SampleID_H1[a]))  


    elif b == "SINGLE LAIV":
       if PigID_H1[a] not in H1_SINGLE_LAIV.keys():
          H1_SINGLE_LAIV[PigID_H1[a]]=[str(SampleID_H1[a])]
       
       else:
          if str(SampleID_H1[a]) not in H1_SINGLE_LAIV[PigID_H1[a]]:
             H1_SINGLE_LAIV[PigID_H1[a]].append(str(SampleID_H1[a])) 
             
    elif b == "NO VAC":
        if PigID_H1[a] not in H1_NO_VAC.keys():
          H1_NO_VAC[PigID_H1[a]]=[str(SampleID_H1[a])]
       
        else:
          if str(SampleID_H1[a]) not in H1_NO_VAC[PigID_H1[a]]:
             H1_NO_VAC[PigID_H1[a]].append(str(SampleID_H1[a]))


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
    
    
  
def get_freq (pig_x, sample_x, mutation_D, variant_dict, freq_dict):

    index_x = variant_dict[pig_x][sample_x].index(mutation_D)
    freq_num = float(freq_dict[pig_x][sample_x][index_x]) 

    return freq_num      

def get_divergence (pig_1, sample_1, pig_2, sample_2, variant_dict, freq_dict):
    
    variant_pool=variant_dict[pig_1][sample_1]+variant_dict[pig_2][sample_2]   
    mutation=list(collections.Counter(variant_pool).keys())       
    count=list(collections.Counter(variant_pool).values())
    mutation_shared=[]        
    for l, i in enumerate(count):
        if int(i) == 2:
           mutation_shared.append(mutation[l])
    total = float(sum(freq_dict[pig_1][sample_1]))+ float(sum(freq_dict[pig_2][sample_2]))       
    if len(mutation_shared) > 0:
       for h in mutation_shared:
           total -= get_freq(pig_1, sample_1, h, variant_dict, freq_dict)
           total -= get_freq(pig_2,sample_2, h, variant_dict, freq_dict)
           total += abs(get_freq(pig_1, sample_1, h, variant_dict, freq_dict) - get_freq(pig_2, sample_2, h, variant_dict, freq_dict))
                  
    return total  

def get_pigID (sampleID, treatment_dict):
    for x in treatment_dict.keys():
        if sampleID in treatment_dict[x]:
           return str(x) 
    
    
def within_treatment (treatment_dict, variant_dict, freq_dict):
    mark=[]
    sample_pool=[]
    distance=[]    
    for x in treatment_dict.keys():
        sample_pool += treatment_dict[x]
        
    for sample_1 in sample_pool:
        for sample_2 in sample_pool:
            if sample_1 != sample_2 and sample_1 + "-" + sample_2 not in mark:
               pig_1 = get_pigID (sample_1, treatment_dict) 
               pig_2 = get_pigID (sample_2, treatment_dict) 
               distance.append(float(get_divergence (pig_1, sample_1, pig_2, sample_2, variant_dict, freq_dict)))
               mark.append(sample_1 +"-"+sample_2)
               mark.append(sample_2+"-"+sample_1)
    return distance           

PRIME_BOOST_H3=within_treatment (H3_PRIME_BOOST, H3_variants, H3_freq)
SINGLE_LAIV_H3=within_treatment (H3_SINGLE_LAIV, H3_variants, H3_freq)
NO_VAC_H3=within_treatment (H3_NO_VAC, H3_variants, H3_freq)

df_H3_PRIME_BOOST = pd.DataFrame({"L1 norm": PRIME_BOOST_H3})
df_H3_PRIME_BOOST = df_H3_PRIME_BOOST.reset_index()
df_H3_PRIME_BOOST["Treatment"]= "PRIME BOOST" 

df_H3_SINGLE_LAIV = pd.DataFrame({"L1 norm": SINGLE_LAIV_H3})
df_H3_SINGLE_LAIV = df_H3_SINGLE_LAIV.reset_index()
df_H3_SINGLE_LAIV["Treatment"]= "SINGLE LAIV" 

df_H3_NO_VAC = pd.DataFrame({"L1 norm": NO_VAC_H3})
df_H3_NO_VAC = df_H3_NO_VAC.reset_index()
df_H3_NO_VAC["Treatment"]= "NO VAC" 


PRIME_BOOST_H1=within_treatment (H1_PRIME_BOOST, H1_variants, H1_freq)
SINGLE_LAIV_H1=within_treatment (H1_SINGLE_LAIV, H1_variants, H1_freq)
NO_VAC_H1=within_treatment (H1_NO_VAC, H1_variants, H1_freq)

df_H1_PRIME_BOOST = pd.DataFrame({"L1 norm": PRIME_BOOST_H1})
df_H1_PRIME_BOOST = df_H1_PRIME_BOOST.reset_index()
df_H1_PRIME_BOOST["Treatment"]= "PRIME BOOST" 

df_H1_SINGLE_LAIV = pd.DataFrame({"L1 norm": SINGLE_LAIV_H1})
df_H1_SINGLE_LAIV = df_H1_SINGLE_LAIV.reset_index()
df_H1_SINGLE_LAIV["Treatment"]= "SINGLE LAIV" 

df_H1_NO_VAC = pd.DataFrame({"L1 norm": NO_VAC_H1})
df_H1_NO_VAC = df_H1_NO_VAC.reset_index()
df_H1_NO_VAC["Treatment"]= "NO VAC" 


def between_treatment (treatment_1_dict, treatment_2_dict, variant_dict, freq_dict):
    mark=[]
    sample_pool_1=[]
    sample_pool_2=[]
    distance=[]  
    for x1 in treatment_1_dict.keys():
        sample_pool_1 += treatment_1_dict[x1]
        
    for x2 in treatment_2_dict.keys():
        sample_pool_2 += treatment_2_dict[x2]    
    
    for sample_1 in sample_pool_1:
        for sample_2 in sample_pool_2:
            if sample_1 + "-" + sample_2 not in mark:
               pig_1 = get_pigID (sample_1, treatment_1_dict) 
               pig_2 = get_pigID (sample_2, treatment_2_dict)
               distance.append(float(get_divergence (pig_1, sample_1, pig_2, sample_2, variant_dict, freq_dict)))
               mark.append(sample_1 +"-"+sample_2)
               mark.append(sample_2+"-"+sample_1)
               
    return distance           
               
PRIME_SINGLE_H3=between_treatment (H3_PRIME_BOOST, H3_SINGLE_LAIV, H3_variants, H3_freq)
SINGLE_NO_H3=between_treatment (H3_SINGLE_LAIV, H3_NO_VAC, H3_variants, H3_freq)
NO_PRIME_H3=between_treatment (H3_NO_VAC, H3_PRIME_BOOST, H3_variants, H3_freq)

df_H3_PRIME_SINGLE = pd.DataFrame({"L1 norm": PRIME_SINGLE_H3})
df_H3_PRIME_SINGLE = df_H3_PRIME_SINGLE.reset_index()
df_H3_PRIME_SINGLE["Treatment"]= "PRIME-SINGLE" 

df_H3_SINGLE_NO = pd.DataFrame({"L1 norm": SINGLE_NO_H3})
df_H3_SINGLE_NO = df_H3_SINGLE_NO.reset_index()
df_H3_SINGLE_NO["Treatment"]= "SINGLE-NO" 

df_H3_NO_PRIME = pd.DataFrame({"L1 norm": NO_PRIME_H3})
df_H3_NO_PRIME = df_H3_NO_PRIME.reset_index()
df_H3_NO_PRIME["Treatment"]= "NO-PRIME" 

df_between_out=pd.concat([df_H3_PRIME_SINGLE,df_H3_SINGLE_NO, df_H3_NO_PRIME, df_H3_PRIME_BOOST,df_H3_SINGLE_LAIV, df_H3_NO_VAC])
           
df_between_out.to_csv(r'/Users/apple/Desktop/H3 treatment L1 norm distance.csv',index=False)       


PRIME_SINGLE_H1=between_treatment (H1_PRIME_BOOST, H1_SINGLE_LAIV, H1_variants, H1_freq)
SINGLE_NO_H1=between_treatment (H1_SINGLE_LAIV, H1_NO_VAC, H1_variants, H1_freq)
NO_PRIME_H1=between_treatment (H1_NO_VAC, H1_PRIME_BOOST, H1_variants, H1_freq)

df_H1_PRIME_SINGLE = pd.DataFrame({"L1 norm": PRIME_SINGLE_H1})
df_H1_PRIME_SINGLE = df_H1_PRIME_SINGLE.reset_index()
df_H1_PRIME_SINGLE["Treatment"]= "PRIME-SINGLE" 

df_H1_SINGLE_NO = pd.DataFrame({"L1 norm": SINGLE_NO_H1})
df_H1_SINGLE_NO = df_H1_SINGLE_NO.reset_index()
df_H1_SINGLE_NO["Treatment"]= "SINGLE-NO" 

df_H1_NO_PRIME = pd.DataFrame({"L1 norm": NO_PRIME_H1})
df_H1_NO_PRIME = df_H1_NO_PRIME.reset_index()
df_H1_NO_PRIME["Treatment"]= "NO-PRIME" 

df_between_out=pd.concat([df_H1_PRIME_SINGLE,df_H1_SINGLE_NO, df_H1_NO_PRIME, df_H1_PRIME_BOOST,df_H1_SINGLE_LAIV, df_H1_NO_VAC])
           
df_between_out.to_csv(r'/Users/apple/Desktop/H1 treatment L1 norm distance.csv',index=False)       


       