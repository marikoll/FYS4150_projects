#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 14:44:27 2018

@author: maritkollstuen
"""

import numpy as np
import numba


def periodic(i, limit, add):
    return (i + limit + add) % limit


@numba.njit(cache = True)
def initial_energy(spins, temp, ordered = False):
    E = 0
    M = 0
    num_spins = len(spins)
    
    for i in range(num_spins):
        for j in range(num_spins):
            if ordered is True: 
                if (temp < 1.5): 
                    spins[i, j] = 1
            left = spins[i-1, j] if i>0 else spins[num_spins - 1, j]
            above = spins[i, j-1] if j>0 else spins[i, num_spins - 1]
            
            E -= spins[i,j]*(left+above)
            M += spins[i, j]
            
    return E, M


@numba.njit(cache=True)
def MC(spins, num_cycles, temperature, ordered = False):#, num_thermalization_steps=0):
    num_spins = len(spins)
    exp_values = np.zeros((num_cycles, 5))
    
    E, M = initial_energy(spins, temperature, ordered)

    for i in range(num_cycles):
        ix = np.random.randint(num_spins)
        iy = np.random.randint(num_spins)

        left = spins[ix - 1, iy] if ix > 0 else spins[num_spins - 1, iy]
        right = spins[ix + 1, iy] if ix < (num_spins - 1) else spins[0, iy]

        above = spins[ix, iy - 1] if iy > 0 else spins[ix, num_spins - 1]
        below = spins[ix, iy + 1] if iy < (num_spins - 1) else spins[ix, 0]

        delta_energy = (2 * spins[ix, iy] * (left + right + above + below))
        if np.random.random() <= np.exp(-delta_energy / temperature):
            spins[ix, iy] *= -1.0
            E += delta_energy
            M += 2*spins[ix, iy]
        
        exp_values[i,0] = E
        exp_values[i,1] = M
        exp_values[i,2] = E**2
        exp_values[i,3] = M**2
        exp_values[i,4] = np.abs(M)

    norm = 1/float(num_cycles)
    
    energy_avg = np.sum(exp_values[:,0])*norm
    magnet_avg = np.sum(exp_values[:,1])*norm
    energy2_avg = np.sum(exp_values[:,2])*norm
    magnet2_avg = np.sum(exp_values[:,3])*norm
    magnet_absavg = np.sum(exp_values[:,4])*norm
    
    energy_var = (energy2_avg - energy_avg**2)/(num_spins**2)#**2*temperature**2)
    magnet_var = (magnet2_avg - magnet_absavg**2)/(num_spins**2)#**2*temperature**2)
    
    energy_avg = energy_avg/num_spins**2
    magnet_avg = magnet_avg/num_spins**2
    C_v = energy_var/temperature**2
    susceptibility = magnet_var/temperature
    abs_magnet = magnet_avg/num_spins

    return energy_avg, magnet_avg, C_v, susceptibility, abs_magnet

    
if __name__ == "__main__": 
    pass
        