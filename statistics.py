#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:40:38 2023

@author: katiagim
"""

import numpy as np
import random
import matplotlib.pyplot as plt
#import collections
import math
import pandas as pd
import statistics
#from matplotlib import animation
import time
#import itertools
import imageio
import os


def mut_position_plot(directory):
    f = pd.read_csv(directory, header=None)
    pp = []
    xtic=[]
    for i in range(len(f.iloc[0])):
        pp.append(np.nansum(f[i]))
    positions=list(range(len(pp)))
    for i in range(len(positions)):
        if pp[i]>=50:
            xtic.append(positions[i])
            
    plt.figure(figsize=(10,10))
    plt.title("Positions in the genome with mutations", size=15)
    plt.xlabel("Position", size=15)
    plt.ylabel("Count", size=15)
    plt.xticks(xtic, rotation=90)
    plt.bar(positions, pp, color='mediumpurple', width=2)
    plt.savefig(directory[:54]+'mutation_positions.jpg', dpi=1200)
    
# def plot_genome(directory):
    
#     f = pd.read_csv(directory, header=None)
    
#     pp = []
#     for i in range(len(f.iloc[0])):
#         pp.append(np.nansum(f[i])) #This counts the mutations in every position in the genome
        
#     for i in range(len(f.iloc[:])):
        
    

directory = "/Users/katiagim/Desktop/MobiVirus/attempt_2/11/genomes/6566.csv"
#mut_position_plot(directory)
n_1 = 0
n_2 = 0
n_3 = 0 
n_4 = 0
n_5 = 0
n_6 = 0

f = pd.read_csv(directory, header=None)
pp = []
for i in range(len(f.iloc[0])):
    pp.append(np.nansum(f[i]))
    
for i in range(1,10):
    if pp[i]>0:
        n_1+=pp[i]
        
for i in range(1, len(pp)):
    if pp[i]>0:
        n_2+=pp[i]
        
for i in range(1, 500):
    if pp[i]>0:
        n_3+=pp[i]    

for i in range(500, 1000):
    if pp[i]>0:
        n_4+=pp[i]    

for i in range(len(f)):
    if f.iloc[i][0]==0 and f.iloc[i][1:10].any()==1:
        n_5+=1
        
for i in range(len(f)):
    if f.iloc[i][0]==1 and f.iloc[i][1:10].any()==1:
        n_6+=1

print("%i Mutations in the first 10 positions" %n_1)
print("%i Mutations in total" %n_2)
print("%i Mutations in the first half of the genome" %n_3)
print("%i Mutations in the second half of the genome" %n_4)
print("%i Mutations from normal spreaders" %n_5)
print("%i Mutations from super spreaders" %n_6 )













