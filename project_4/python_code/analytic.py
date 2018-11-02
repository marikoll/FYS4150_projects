#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 12:07:38 2018

@author: maritkollstuen
"""
import numpy as np
import math


temp_init   = 1.8
temp_end    = 2.6
temp_step   = 0.1


temps = 1 #np.arange(temp_init,temp_end+temp_step/2,temp_step,float)
k_b = 1


configs = [-8, 0, 0, 0, 0,0, 0, 0, 0, 8, 8, 0, 0, 0, 0, -8]


#partition = np.zeros((len(temps), len(configs)), np.float64)
#probability = np.zeros((len(temps), len(configs)), np.float64)
#partition_func = np.zeros(len(temps), np.float64)


#for i in range(len(temps)): 
#    for j in range(len(configs)): 
#        partition[i, j] = math.exp(-(1/temps[i]*k_b)*float(configs[j]))
#    partition_func[i] = sum(partition[i,:])
#    for j in range(len(configs)):    
#        probability[i, j] = math.exp(-(1/temps[i]*k_b)*float(configs[j]))/partition_func[i]
#   

partition = np.zeros(len(configs), np.float64)
probability = np.zeros(len(configs), np.float64)

for j in range(len(configs)): 
    partition[j] = math.exp(-(1/temps*k_b)*float(configs[j]))
z = sum(partition)
for j in range(len(configs)):    
    probability[j] = math.exp(-(1/temps*k_b)*float(configs[j]))/z
    
mean_E = (1/z)*sum([i**2*math.exp(-i) for i in configs])
c_v  = (1/z)*sum([i**2*math.exp(-i) for i in configs]) - ((1/z)*sum([i*math.exp(-i) for i in configs]))**2
mean_M = np.zeros(len(configs), np.float64)
suscept = np.zeros(len(configs), np.float64)



