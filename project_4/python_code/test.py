#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 09:55:35 2018

@author: maritkollstuen
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from numba import jit 
import time

def periodic(i, limit, add):
    return (i + limit + add) % limit

def initialize(spin_matrix, temp, E, M):
    size = len(spin_matrix)
    for y in range(size):
        for x in range(size):
            if temp <1.5:
                spin_matrix[x, y] = 1
            M +=spin_matrix[x, y]
    for y in range(size):
        for x in range(size):
            E -= spin_matrix.item(x,y)*\
                 (spin_matrix.item(periodic(x,size,-1),y) + \
                  spin_matrix.item(x,periodic(y,size,1)))
    
    return E, M

@jit#(nopython = True)    
def Monte_Carlo(temp, size, trials):
    spin_matrix = np.ones((size, size), np.int8)
    E_avg = M_avg = 0
    E2_avg = M2_avg = 0
    # Possible energy changes:
    w = np.zeros(17,np.float64)
    for de in range(-8,9,4): #include +8
        w[de+8] = math.exp(-de/temp)
    
    # Initial energy: 
    E, M = initialize(spin_matrix, temp, E=0, M=0)
    
    
    for i in range(trials):
        for j in range(size**2):
            x = int(np.random.random()*size)
            y = int(np.random.random()*size)
            deltaE = 2*spin_matrix.item(x,y)*\
                     (spin_matrix.item(periodic(x,size,-1), y) +\
                      spin_matrix.item(periodic(x,size,1),  y) +\
                      spin_matrix.item(x, periodic(y,size,-1)) +\
                      spin_matrix.item(x, periodic(y,size,1)))
            if np.random.random() <= w[deltaE+8]:
                #Accept!
                spin_matrix[x,y] *= -1
                M += 2*spin_matrix[x, y]
                E += deltaE
            
        #Update expectation values
        E_avg   += E
        E2_avg  += E**2
        M_avg   += M
        M2_avg  += M**2
        

    E_avg   /= float(trials)
    E2_avg  /= float(trials)
    M_avg   /= float(trials) 
    M2_avg  /= float(trials)
    #Calculate variance and normalize to per-point and temp
    E_variance  = (E2_avg - E_avg**2)/float(size**2*temp**2)
    M_variance  = (M2_avg - M_avg**2)/float(size**2*temp**2) #susceptibility
    #Normalize returned averages to per-point
    E_avg /= float(size**2)
    M_avg /= float(size**2)

    return (E_avg, E_variance, M_avg, M_variance)
    


if __name__ == "__main__": 
    size        = 2
    trials      = 1000
    temp_init   = 1.8
    temp_end    = 2.6
    temp_step   = 0.1
    
    
    temps = np.arange(temp_init,temp_end+temp_step/2,temp_step,float)
    temps = [float(1)]
    Dim = np.size(temps)
    energy = np.zeros(Dim)
    heatcapacity = np.zeros(Dim) 
    temperature = np.zeros(Dim)
    magnetization = np.zeros(Dim)
    susceptibility = np.zeros(Dim)
    start = time.time()
    for i in range(Dim):
        (E_avg, E_variance, M_avg, M_variance) = Monte_Carlo(temps[i],size,trials)
        temperature[i] = temps[i]
        energy[i] = E_avg
        heatcapacity[i] = E_variance
        magnetization[i] = M_avg
        susceptibility[i] = M_variance
    stop = time.time()
    print('Monte Carlo: {:.3f} s'.format(stop - start))
    
#    plt.figure(1)
#    plt.subplot(211)
##    plt.axis([1.8,2.6,-2.0, -1.0])
#    plt.xlabel(r'Temperature $J/(k_B)$')
#    plt.ylabel(r'Average energy per spin  $E/N$')
#    plt.plot(temperature, energy, 'b-')
#    plt.subplot(212)
##    plt.axis([1.8,2.6, 0.0, 2.0])
#    plt.plot(temperature, heatcapacity, 'r-')
#    plt.xlabel(r'Temperature $J/(k_B)$')
#    plt.ylabel(r'Heat capacity per spin  $C_V/N$')
##    plt.savefig('energycv.pdf')
#    plt.show()
#    
#    plt.figure(2)
#    plt.subplot(211)
##    plt.axis([1.8,2.6,-2.0, -1.0])
#    plt.xlabel(r'Temperature $J/(k_B)$')
#    plt.ylabel(r'Average energy per spin  $E/N$')
#    plt.plot(temperature, magnetization, 'b-')
#    plt.subplot(212)
##    plt.axis([1.8,2.6, 0.0, 2.0])
#    plt.plot(temperature, susceptibility, 'r-')
#    plt.xlabel(r'Temperature $J/(k_B)$')
#    plt.ylabel(r'Heat capacity per spin  $C_V/N$')
##    plt.savefig('energycv.pdf')
#    plt.show()